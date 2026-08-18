from pipeline.clients.data_loder import DataLoader


def load_client_data(partition_id: int) -> DataLoader:
    """Load client configuration for the given partition ID."""
    config_file_path = f"configs/center_{partition_id}/config.yaml"
    return DataLoader(config_file_path)
