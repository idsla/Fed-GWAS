# server/main_server.py

from flwr.common import Context
from flwr.server import ServerApp, ServerAppComponents, ServerConfig
import logging

# Use the strict strategy implementation (Scalar-only config, safer with Flower)
from pipeline.server.strategy_strict import FederatedGWASStrategy
from pipeline.utils.monitoring_config import resolve_monitoring_settings
from pipeline.utils.retention_config import resolve_retention_settings

logger = logging.getLogger(__name__)


def _resolve_server_config_file(config_path):
    """Resolve server config for both CLI flat layout and legacy nested layout.

    New CLI projects write `configs/server.yaml`. Existing experiments still
    often provide `configs/server/config.yaml`. Resolving both here keeps
    `fedgwas-sim run` compatible without changing older experiment folders.
    """
    from pathlib import Path

    config_path_obj = Path(config_path)
    if config_path_obj.suffix in (".yaml", ".yml"):
        return config_path_obj
    if config_path_obj.name == "configs":
        candidates = [
            config_path_obj / "server.yaml",
            config_path_obj / "server" / "config.yaml",
        ]
    else:
        candidates = [
            config_path_obj / "configs" / "server.yaml",
            config_path_obj / "configs" / "server" / "config.yaml",
        ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def server_fn(context: Context):

    num_rounds = context.run_config["num-server-rounds"]
    import yaml
    from pathlib import Path

    def _read_int_run_config(keys, default):
        for key in keys:
            if key not in context.run_config:
                continue
            raw_value = context.run_config[key]
            try:
                return int(raw_value)
            except (TypeError, ValueError):
                logger.warning(
                    "[Server App] Invalid run_config value for %s=%r; using default=%s",
                    key,
                    raw_value,
                    default,
                )
                return default
        return default

    def _read_int_from_mapping(mapping, keys, default):
        value = mapping
        for key in keys:
            if not isinstance(value, dict) or key not in value:
                return default
            value = value[key]
        try:
            return int(value)
        except (TypeError, ValueError):
            return default

    # Server output directory must be declared in server config file
    output_dir = None
    
    # Priority 1: Check for explicit server output_dir in run_config (for programmatic override)
    if "server_output_dir" in context.run_config:
        output_dir = context.run_config["server_output_dir"]
        logger.info("[Server App] Using explicit server output directory from run_config: %s", output_dir)
    else:
        # Priority 2: Load from server config file: configs/server/config.yaml
        config_path = context.run_config.get("config_path")
        if not config_path:
            raise ValueError("config_path not found in run_config. Cannot locate server config file.")
        
        server_config_file = _resolve_server_config_file(config_path)
        
        if not server_config_file.exists():
            raise ValueError(f"Server config file not found: {server_config_file}. Server output directory must be declared in configs/server/config.yaml")
        
        try:
            with open(server_config_file, 'r') as f:
                server_config = yaml.safe_load(f)
            # Server config must have output.log_dir or output.intermediate_dir
            if "output" not in server_config:
                raise ValueError(f"Server config file {server_config_file} must have 'output' section")
            
            if "log_dir" in server_config["output"]:
                # Extract server directory from log_dir (e.g., .../results/server/logs -> .../results/server)
                log_dir = Path(server_config["output"]["log_dir"])
                if log_dir.name == "logs" and log_dir.parent.name == "server":
                    output_dir = str(log_dir.parent)
                else:
                    # Try to construct server directory from log_dir parent
                    output_dir = str(log_dir.parent / "server") if log_dir.parent.name != "server" else str(log_dir.parent)
                logger.info("[Server App] Using server output directory from server config: %s", output_dir)
            elif "intermediate_dir" in server_config["output"]:
                # Extract server directory from intermediate_dir
                interm_dir = Path(server_config["output"]["intermediate_dir"])
                if interm_dir.parent.name == "server":
                    output_dir = str(interm_dir.parent)
                else:
                    output_dir = str(interm_dir.parent / "server")
                logger.info("[Server App] Using server output directory from server config: %s", output_dir)
            else:
                raise ValueError(f"Server config file {server_config_file} must have output.log_dir or output.intermediate_dir")
        except ValueError:
            raise
        except Exception as e:
            raise ValueError(f"Failed to load server config file {server_config_file}: {e}")
    
    if output_dir is None:
        raise ValueError("Server output directory not configured. Must be set via server_output_dir in run_config or declared in configs/server/config.yaml")

    # Resolve defaults from experiment-level config.yaml when available
    default_chunk_size = 1000
    default_lr_pad_to = 0
    config_path_value = context.run_config.get("config_path")
    if config_path_value:
        config_path_obj = Path(config_path_value)
        if config_path_obj.suffix in (".yaml", ".yml"):
            experiment_config_file = config_path_obj
        elif config_path_obj.name == "configs":
            experiment_config_file = config_path_obj.parent / "config.yaml"
        else:
            experiment_config_file = config_path_obj / "config.yaml"

        if experiment_config_file.exists():
            try:
                with open(experiment_config_file, "r") as f:
                    experiment_config = yaml.safe_load(f) or {}
                default_chunk_size = _read_int_from_mapping(
                    experiment_config, ("server", "chunk_size"), default_chunk_size
                )
                default_lr_pad_to = _read_int_from_mapping(
                    experiment_config, ("server", "lr_pad_to"), default_lr_pad_to
                )
            except Exception as e:
                logger.warning(
                    "[Server App] Failed to read experiment defaults from %s: %s",
                    experiment_config_file,
                    e,
                )

    chunk_size = _read_int_run_config(
        ("chunk_size", "chunk-size", "server_chunk_size", "server-chunk-size"),
        default=default_chunk_size,
    )
    lr_pad_to = _read_int_run_config(("lr_pad_to", "lr-pad-to"), default=default_lr_pad_to)
    if chunk_size < 1:
        logger.warning("[Server App] chunk_size=%s is invalid; forcing chunk_size=1", chunk_size)
        chunk_size = 1
    if lr_pad_to < 0:
        logger.warning("[Server App] lr_pad_to=%s is invalid; forcing lr_pad_to=0", lr_pad_to)
        lr_pad_to = 0

    monitoring_settings = resolve_monitoring_settings(
        config_path=str(config_path_value) if config_path_value else None
    )
    retention_settings = resolve_retention_settings(
        config_path=str(config_path_value) if config_path_value else None
    )

    strategy = FederatedGWASStrategy(
        output_dir=output_dir,
        chunk_size=chunk_size,
        lr_pad_to=lr_pad_to,
        monitoring_settings=monitoring_settings,
        retention_settings=retention_settings,
    )

    config = ServerConfig(num_rounds=num_rounds)
    return ServerAppComponents(config=config, strategy=strategy)


app = ServerApp(server_fn=server_fn)
