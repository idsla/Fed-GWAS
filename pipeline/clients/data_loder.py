# client/data_loader.py

import os
import yaml
import subprocess
import logging

class DataLoader:
    
    def __init__(
        self, 
        config_file="config.yaml"
    ):
        
        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)
        
        self.input_path = self.config["input_data"]["path"]
        self.input_type = self.config["input_data"].get("type", "bed")
        self.intermediate_dir = self.config["output"].get("intermediate_dir", "intermediate")
        self.log_dir = self.config["output"].get("log_dir", "logs")
        self.thresholds = self.config.get("thresholds", {})
        self.parameters = self.config.get("parameters", {})
        
        # Participation flags for each pipeline stage
        self.participation = self.config.get("participation", {})

        self.flower_config = self.config.get("flower", {})

        os.makedirs(self.intermediate_dir, exist_ok=True)
        os.makedirs(self.log_dir, exist_ok=True)

    def transform_data(self):
        """
        Transform the input data into PLINK binary format if needed.
        If the input type is 'bed', assume it's already in the correct format.
        If the input type is 'vcf', convert it using PLINK.
        Returns the dataset prefix (without extension).
        """
        if self.input_type.lower() == "bed":
            return self.input_path
        elif self.input_type.lower() == "vcf":
            out_prefix = f"{self.intermediate_dir}/converted_data"
            cmd = [
                "plink",
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

