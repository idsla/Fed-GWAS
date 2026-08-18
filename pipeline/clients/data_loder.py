# client/data_loader.py

import logging
import os
import subprocess
from pathlib import Path
from typing import Dict

import yaml
from flwr.common import ConfigRecord


def _repo_root_from_config(config_file: str) -> Path:
    """Find repo root (pyproject.toml) from an absolute client config path."""
    start = Path(config_file).resolve()
    for parent in [start.parent, *start.parents]:
        if (parent / "pyproject.toml").is_file():
            return parent
    return Path.cwd()


def _abs_repo_path(path: str, repo_root: Path) -> str:
    p = Path(path)
    if p.is_absolute():
        return str(p.resolve())
    return str((repo_root / p).resolve())


class DataLoader:
    
    def __init__(
        self, 
        config_file: str
    ):
        config_file = str(Path(config_file).resolve())
        self.config_file = config_file
        self.repo_root = _repo_root_from_config(config_file)
        
        with open(config_file, "r", encoding="utf-8") as f:
            self.config = yaml.safe_load(f)
        
        self.input_path = _abs_repo_path(self.config["input_data"]["path"], self.repo_root)
        self.input_type = self.config["input_data"].get("type", "bed")
        self.intermediate_dir = _abs_repo_path(
            self.config["output"].get("intermediate_dir", "intermediate"), self.repo_root
        )
        self.log_dir = _abs_repo_path(
            self.config["output"].get("log_dir", "logs"), self.repo_root
        )
        self.thresholds = self.config.get("thresholds", {})
        self.parameters = self.config.get("parameters", {})
        
        # create intermediate and log directories
        os.makedirs(self.intermediate_dir, exist_ok=True)
        os.makedirs(self.log_dir, exist_ok=True)
        
        # Participation flags for each pipeline stage
        self.participation = self.config.get("participation", {})
        self.flower_config = self.config.get("flower", {})
        self.monitoring = self.config.get("monitoring", {})
        
        # prefix for plink files
        self.plink_prefix = None

    def transform_data(self):
        """
        Transform the input data into PLINK binary format if needed.
        If the input type is 'bed', assume it's already in the correct format.
        If the input type is 'vcf', convert it using PLINK.
        Returns the dataset prefix (without extension).
        """
        if self.input_type.lower() == "bed":
            self.plink_prefix = self.input_path
            bed = f"{self.plink_prefix}.bed"
            if not os.path.isfile(bed):
                raise FileNotFoundError(
                    f"PLINK bed not found: {bed} "
                    f"(repo_root={self.repo_root}, config={self.config_file})"
                )
            return self.plink_prefix
        elif self.input_type.lower() == "vcf":
            # Ensure intermediate directory exists for VCF conversion
            os.makedirs(self.intermediate_dir, exist_ok=True)
            out_prefix = os.path.join(self.intermediate_dir, "converted_data")
            cmd = [
                "./bin/plink",
                "--vcf", self.input_path,
                "--make-bed",
                "--out", out_prefix
            ]
            try:
                subprocess.run(cmd, check=True)
                self.plink_prefix = out_prefix
                logging.info(f"Converted VCF {self.input_path} to PLINK binary {out_prefix}.")
                return self.plink_prefix
            except subprocess.CalledProcessError as e:
                logging.error(f"Error converting VCF: {e}")
                raise
        else:
            raise ValueError(f"Unsupported input type: {self.input_type}")

    def get_thresholds(self):
        return self.thresholds

    def get_parameters(self):
        return self.parameters

    def get_flower_config(self):
        return self.flower_config

    def get_participation(self):
        """
        Return the participation flags dict indicating which stages the client participates in.
        """
        return self.participation

    def get_monitoring(self):
        """Return monitoring flags (may be empty; merged from experiment config at runtime)."""
        return self.monitoring

    
    def to_config_records(self):
        """
        Convert all DataLoader attributes to a dictionary suitable for ConfigRecords.
            
        Returns:
            Dict[str, ConfigRecord]: Dictionary containing all DataLoader configuration
        """
        config_records = {}
        
        config_records['client_paths'] =  ConfigRecord({
            'input_path': self.input_path,
            'input_type': self.input_type,
            'plink_prefix': self.plink_prefix,
            "intermediate_dir": self.intermediate_dir,
            "log_dir": self.log_dir
        })
        
        config_records['thresholds'] = ConfigRecord(self.thresholds)
        
        config_records['parameters'] = ConfigRecord(self.parameters)
        
        config_records['participation'] = ConfigRecord(self.participation)
        
        config_records['flower_config'] = ConfigRecord(self.flower_config)

        if self.monitoring:
            config_records['monitoring'] = ConfigRecord(self.monitoring)
            
        return config_records
    
    @staticmethod
    def from_config_records(config_records: Dict[str, ConfigRecord]):
        """
        Create a DataLoader-like object from stored configuration.
        
        Args:
            config_dict (Dict[str, ConfigRecord]): Dictionary containing stored configuration
            
        Returns:
            DataLoader: DataLoader instance with restored configuration
        """
        # Create a dummy DataLoader without reading config file
        loader = DataLoader.__new__(DataLoader)
        
        # Restore all attributes from config_dict
        loader.input_path = config_records['client_paths']["input_path"]
        loader.input_type = config_records['client_paths']["input_type"]
        loader.plink_prefix = config_records['client_paths']["plink_prefix"]
        loader.intermediate_dir = config_records['client_paths']["intermediate_dir"]
        loader.log_dir = config_records['client_paths']["log_dir"]
        loader.thresholds = config_records['thresholds']
        loader.parameters = config_records['parameters']
        loader.participation = config_records['participation']
        loader.flower_config = config_records['flower_config']
        loader.monitoring = dict(config_records['monitoring']) if 'monitoring' in config_records else {}
        
        return loader
    
    def __str__(self):
        data_loader_desc = "DataLoader Object:\n"
        data_loader_desc += f"  Input Path: {self.input_path} Input Type: {self.input_type}\n"
        data_loader_desc += f"  Intermediate Dir: {self.intermediate_dir}\n"
        data_loader_desc += f"  Log Dir: {self.log_dir}\n"
        data_loader_desc += f"  Thresholds: {self.thresholds}\n"
        data_loader_desc += f"  Parameters: {self.parameters}\n"
        data_loader_desc += f"  Participation: {self.participation}\n"
        data_loader_desc += f"  Flower Config: {self.flower_config}\n"
        return data_loader_desc
