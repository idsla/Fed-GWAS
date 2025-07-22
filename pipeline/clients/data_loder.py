# client/data_loader.py

import os
import yaml
import subprocess
import logging
from flwr.common import ConfigRecord, RecordDict
from typing import Dict

class DataLoader:
    
    def __init__(
        self, 
        config_file: str
    ):
        
        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)
        
        self.input_path = self.config["input_data"]["path"]
        self.input_type = self.config["input_data"].get("type", "bed")
        self.intermediate_dir = self.config["output"].get("intermediate_dir", "intermediate")
        self.log_dir = self.config["output"].get("log_dir", "logs")
        self.thresholds = self.config.get("thresholds", {})
        self.parameters = self.config.get("parameters", {})
        
        # create intermediate and log directories
        os.makedirs(self.intermediate_dir, exist_ok=True)
        os.makedirs(self.log_dir, exist_ok=True)
        
        # Participation flags for each pipeline stage
        self.participation = self.config.get("participation", {})
        self.flower_config = self.config.get("flower", {})
        
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
            return self.plink_prefix
        elif self.input_type.lower() == "vcf":
            # TODO: check (this branch does not have return)
            # Ensure intermediate directory exists for VCF conversion
            os.makedirs(self.intermediate_dir, exist_ok=True)
            out_prefix = f"{self.intermediate_dir}/converted_data"
            cmd = [
                "./bin/plink",
                "--vcf", self.input_path,
                "--make-bed",
                "--out", out_prefix
            ]
            try:
                subprocess.run(cmd, check=True)
                logging.info(f"Converted VCF {self.input_path} to PLINK binary {out_prefix}.")
                return out_prefix
            except subprocess.CalledProcessError as e:
                logging.error(f"Error converting VCF: {e}")
                raise e
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

