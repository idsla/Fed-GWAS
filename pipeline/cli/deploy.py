"""Typer command handlers for FedGWAS deployment mode.

The deployment CLI wraps the Flower deployment runtime used by the legacy
cluster scripts: SuperLink on the server, one SuperNode per center, and
`flwr run` submitted through the SuperLink Exec API.
"""

from __future__ import annotations

import re
import shutil
import socket
import subprocess
import sys
from enum import Enum
from pathlib import Path
from typing import Annotated, Optional

import typer


class ScaleChoice(str, Enum):
    """Repository experiment scale shortcuts accepted by `server run`."""

    tiny = "tiny"
    small = "small"
    medium = "medium"


ROOT_EXAMPLES = """Examples:\n
  fedgwas-deploy server start --host 0.0.0.0 --daemon\n
  fedgwas-deploy client start --server 192.168.1.88 --center-id 1 --config configs/center_1.yaml --daemon\n
  fedgwas-deploy server run --server 192.168.1.88 --rounds 20 --scale tiny\n
  fedgwas-deploy server status\n
  fedgwas-deploy client stop\n
"""

app = typer.Typer(
    help="Run FedGWAS in Flower deployment mode.",
    epilog=ROOT_EXAMPLES,
    no_args_is_help=True,
)
server_app = typer.Typer(help="Manage the deployment server node.", no_args_is_help=True)
client_app = typer.Typer(help="Manage a deployment client center node.", no_args_is_help=True)

app.add_typer(server_app, name="server")
app.add_typer(client_app, name="client")


SCALE_CONFIG_PATHS = {
    ScaleChoice.tiny: Path("experiments/correctness/tiny_even/configs"),
    ScaleChoice.small: Path("experiments/performance/small_even/configs"),
    ScaleChoice.medium: Path("experiments/performance/medium_even/configs"),
}


def _run_or_daemon(
    cmd: list[str],
    *,
    daemon: bool,
    log_file: Optional[Path],
    cwd: Path,
) -> int:
    """Run a command in the foreground or start it in the background."""
    if not daemon:
        result = subprocess.run(cmd, cwd=cwd, check=False)
        return int(result.returncode)

    resolved_log = log_file or Path("/tmp") / f"{Path(cmd[0]).name}.log"
    resolved_log.parent.mkdir(parents=True, exist_ok=True)
    log_handle = resolved_log.open("ab")
    process = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        start_new_session=(sys.platform != "win32"),
    )
    typer.echo(f"Started {' '.join(cmd[:1])} pid={process.pid} log={resolved_log}")
    return 0


def _format_host_port(host: str, port: int) -> str:
    """Return Flower-compatible host:port text."""
    if ":" in host and not host.startswith("["):
        return f"[{host}]:{port}"
    return f"{host}:{port}"


def _default_client_port(center_id: int) -> int:
    """Return the legacy cluster script default ClientAppIo port."""
    return 9093 + center_id


def _resolve_config_path(config: Path, cwd: Path) -> Path:
    """Resolve a center config path relative to the current working directory."""
    candidate = config if config.is_absolute() else cwd / config
    return candidate.resolve()


def _scale_config_path(scale: ScaleChoice, cwd: Path) -> Path:
    """Resolve a repository scale shortcut to a config directory."""
    return (cwd / SCALE_CONFIG_PATHS[scale]).resolve()


def _display_path(path: Path, cwd: Path) -> str:
    """Return a stable project-relative path when possible."""
    try:
        return path.resolve().relative_to(cwd.resolve()).as_posix()
    except ValueError:
        return str(path)


def _patch_federation_addresses(pyproject: Path, address: str) -> str:
    """Patch deployment federation addresses and return the original file text."""
    original = pyproject.read_text(encoding="utf-8")
    text = original
    for section in ("local-deployment", "cluster-deployment"):
        pattern = (
            rf'(\[tool\.flwr\.federations\.{re.escape(section)}\][^\[]*?'
            rf'address = ")[^"]+(")'
        )
        text, count = re.subn(
            pattern,
            lambda match, value=address: f"{match.group(1)}{value}{match.group(2)}",
            text,
            count=1,
            flags=re.DOTALL,
        )
        if count == 0:
            raise typer.BadParameter(
                f"Could not patch address in [tool.flwr.federations.{section}] in {pyproject}"
            )
    pyproject.write_text(text, encoding="utf-8")
    return original


def _check_command(name: str) -> tuple[bool, str]:
    found = shutil.which(name)
    return found is not None, found or "not found on PATH"


def _check_tcp(host: str, port: int, timeout: float = 2.0) -> tuple[bool, str]:
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True, f"{host}:{port} reachable"
    except OSError as exc:
        return False, f"{host}:{port} unreachable ({exc})"


def _print_checks(rows: list[tuple[str, bool, str]]) -> int:
    failures = 0
    for label, ok, detail in rows:
        status = "PASS" if ok else "FAIL"
        typer.echo(f"{status} {label}: {detail}")
        if not ok:
            failures += 1
    return failures


SERVER_START_EXAMPLES = """Examples:\n
  fedgwas-deploy server start --host 0.0.0.0\n
  fedgwas-deploy server start --host 0.0.0.0 --daemon --log-file /tmp/superlink.log\n
"""


@server_app.command("start", help="Start Flower SuperLink on the server node.", epilog=SERVER_START_EXAMPLES)
def server_start(
    host: Annotated[str, typer.Option("--host", help="Address SuperLink should bind to.")] = "0.0.0.0",
    fleet_port: Annotated[int, typer.Option("--fleet-port", help="Fleet API port for SuperNodes.")] = 9092,
    exec_port: Annotated[int, typer.Option("--exec-port", help="Exec API port for flwr run.")] = 9093,
    insecure: Annotated[bool, typer.Option("--insecure/--secure", help="Use Flower insecure transport.")] = True,
    daemon: Annotated[bool, typer.Option("--daemon", help="Start in background and write a log file.")] = False,
    log_file: Annotated[Optional[Path], typer.Option("--log-file", help="Daemon log file path.")] = None,
) -> None:
    cmd = ["flower-superlink"]
    if insecure:
        cmd.append("--insecure")
    cmd.extend(
        [
            "--fleet-api-address",
            _format_host_port(host, fleet_port),
            "--exec-api-address",
            _format_host_port(host, exec_port),
        ]
    )
    raise typer.Exit(_run_or_daemon(cmd, daemon=daemon, log_file=log_file, cwd=Path.cwd()))


SERVER_RUN_EXAMPLES = """Examples:\n
  fedgwas-deploy server run --server 192.168.1.88 --rounds 20 --scale tiny\n
  fedgwas-deploy server run --server 192.168.1.88 --rounds 80 --config-path experiments/performance/small_even/configs\n
  fedgwas-deploy server run --federation local-deployment --rounds 20 --config-path configs --no-stream\n
"""


@server_app.command("run", help="Submit the FedGWAS app through Flower Exec API.", epilog=SERVER_RUN_EXAMPLES)
def server_run(
    server: Annotated[
        Optional[str],
        typer.Option("--server", help="SuperLink server host used to patch pyproject federation addresses."),
    ] = None,
    exec_port: Annotated[int, typer.Option("--exec-port", help="SuperLink Exec API port.")] = 9093,
    federation: Annotated[str, typer.Option("--federation", help="Flower federation name.")] = "local-deployment",
    rounds: Annotated[int, typer.Option("--rounds", "-r", help="Number of Flower server rounds.")] = 100,
    scale: Annotated[ScaleChoice, typer.Option("--scale", help="Repository config shortcut.")] = ScaleChoice.tiny,
    config_path: Annotated[Optional[Path], typer.Option("--config-path", help="Config directory to pass to FedGWAS.")] = None,
    stream: Annotated[bool, typer.Option("--stream/--no-stream", help="Stream Flower run logs.")] = True,
) -> None:
    cwd = Path.cwd()
    resolved_config = (cwd / config_path).resolve() if config_path else _scale_config_path(scale, cwd)
    if not resolved_config.is_dir():
        raise typer.BadParameter(f"Config path not found: {_display_path(resolved_config, cwd)}")

    config_display = _display_path(resolved_config, cwd)
    cmd = ["flwr", "run", ".", federation]
    if stream:
        cmd.append("--stream")
    cmd.extend(
        [
            "--run-config",
            f'simulation=false num-server-rounds={rounds} config_path="{config_display}"',
        ]
    )

    pyproject = cwd / "pyproject.toml"
    original_text: str | None = None
    try:
        if server:
            if not pyproject.is_file():
                raise typer.BadParameter("pyproject.toml not found in the current directory")
            original_text = _patch_federation_addresses(pyproject, _format_host_port(server, exec_port))
        result = subprocess.run(cmd, cwd=cwd, check=False)
        raise typer.Exit(result.returncode)
    finally:
        if original_text is not None:
            pyproject.write_text(original_text, encoding="utf-8")


@server_app.command("check", help="Check server-side deployment prerequisites.")
def server_check(
    host: Annotated[str, typer.Option("--host", help="Host to probe for listening ports.")] = "127.0.0.1",
    fleet_port: Annotated[int, typer.Option("--fleet-port", help="Fleet API port.")] = 9092,
    exec_port: Annotated[int, typer.Option("--exec-port", help="Exec API port.")] = 9093,
    network: Annotated[bool, typer.Option("--network/--no-network", help="Probe SuperLink ports.")] = False,
) -> None:
    rows = [
        ("flower-superlink", *_check_command("flower-superlink")),
        ("flwr", *_check_command("flwr")),
        ("pyproject.toml", (Path.cwd() / "pyproject.toml").is_file(), str(Path.cwd() / "pyproject.toml")),
    ]
    if network:
        rows.append(("Fleet API", *_check_tcp(host, fleet_port)))
        rows.append(("Exec API", *_check_tcp(host, exec_port)))
    raise typer.Exit(1 if _print_checks(rows) else 0)


@server_app.command("status", help="Show local Flower deployment processes and ports.")
def server_status() -> None:
    cmd = ["cmd", "/c", "netstat -ano | findstr \"9092 9093 9094 9095\""] if sys.platform == "win32" else [
        "bash",
        "-lc",
        "ps aux | grep -E 'flower-superlink|flower-supernode|flwr run' | grep -v grep; "
        "(netstat -tuln 2>/dev/null || ss -tuln 2>/dev/null) | grep -E '9092|9093|9094|9095' || true",
    ]
    raise typer.Exit(subprocess.run(cmd, check=False).returncode)


@server_app.command("stop", help="Stop local Flower deployment processes on this node.")
def server_stop() -> None:
    if sys.platform == "win32":
        cmd = ["powershell", "-NoProfile", "-Command", "Get-Process | ? {$_.ProcessName -match 'flower|flwr'} | Stop-Process"]
    else:
        cmd = ["bash", "-lc", "pkill -f 'flower-superlink' || true; pkill -f 'flower-supernode' || true; pkill -f 'flwr run' || true"]
    raise typer.Exit(subprocess.run(cmd, check=False).returncode)


CLIENT_START_EXAMPLES = """Examples:\n
  fedgwas-deploy client start --server 192.168.1.88 --center-id 1 --config configs/center_1.yaml\n
  fedgwas-deploy client start --server 192.168.1.88 --center-id 2 --config configs/center_2.yaml --daemon\n
"""


@client_app.command("start", help="Start Flower SuperNode for one center.", epilog=CLIENT_START_EXAMPLES)
def client_start(
    server: Annotated[str, typer.Option("--server", help="SuperLink server host.")],
    center_id: Annotated[int, typer.Option("--center-id", help="One-based center id.")],
    config: Annotated[Path, typer.Option("--config", help="Center config YAML file.")],
    num_centers: Annotated[int, typer.Option("--num-centers", help="Total number of centers/SuperNodes.")] = 2,
    port: Annotated[Optional[int], typer.Option("--port", help="ClientAppIo port. Defaults to 9093 + center id.")] = None,
    host: Annotated[str, typer.Option("--host", help="ClientAppIo bind address.")] = "0.0.0.0",
    fleet_port: Annotated[int, typer.Option("--fleet-port", help="SuperLink Fleet API port.")] = 9092,
    insecure: Annotated[bool, typer.Option("--insecure/--secure", help="Use Flower insecure transport.")] = True,
    daemon: Annotated[bool, typer.Option("--daemon", help="Start in background and write a log file.")] = False,
    log_file: Annotated[Optional[Path], typer.Option("--log-file", help="Daemon log file path.")] = None,
) -> None:
    if center_id < 1:
        raise typer.BadParameter("--center-id must be >= 1")
    if center_id > num_centers:
        raise typer.BadParameter("--center-id cannot exceed --num-centers")

    cwd = Path.cwd()
    resolved_config = _resolve_config_path(config, cwd)
    if not resolved_config.is_file():
        raise typer.BadParameter(f"Center config file not found: {resolved_config}")

    client_port = port or _default_client_port(center_id)
    partition_id = center_id - 1
    cmd = ["flower-supernode"]
    if insecure:
        cmd.append("--insecure")
    cmd.extend(
        [
            "--superlink",
            _format_host_port(server, fleet_port),
            "--clientappio-api-address",
            _format_host_port(host, client_port),
            "--node-config",
            f'partition-id={partition_id} num-partitions={num_centers} config-file="{resolved_config}"',
        ]
    )
    raise typer.Exit(_run_or_daemon(cmd, daemon=daemon, log_file=log_file, cwd=cwd))


@client_app.command("check", help="Check client-side deployment prerequisites.")
def client_check(
    server: Annotated[str, typer.Option("--server", help="SuperLink server host.")],
    center_id: Annotated[int, typer.Option("--center-id", help="One-based center id.")],
    config: Annotated[Path, typer.Option("--config", help="Center config YAML file.")],
    fleet_port: Annotated[int, typer.Option("--fleet-port", help="SuperLink Fleet API port.")] = 9092,
    network: Annotated[bool, typer.Option("--network/--no-network", help="Probe server Fleet API.")] = False,
) -> None:
    cwd = Path.cwd()
    resolved_config = _resolve_config_path(config, cwd)
    rows = [
        ("flower-supernode", *_check_command("flower-supernode")),
        ("center id", center_id >= 1, str(center_id)),
        ("config file", resolved_config.is_file(), str(resolved_config)),
        ("plink", *_check_command("plink")),
    ]
    if network:
        rows.append(("Fleet API", *_check_tcp(server, fleet_port)))
    raise typer.Exit(1 if _print_checks(rows) else 0)


@client_app.command("status", help="Show local Flower deployment processes and ports.")
def client_status() -> None:
    server_status()


@client_app.command("stop", help="Stop local Flower deployment processes on this node.")
def client_stop() -> None:
    server_stop()


def main() -> None:
    """Run the deployment Typer app."""
    app()
