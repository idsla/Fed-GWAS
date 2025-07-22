
from pipeline.clients.data_loder import DataLoader


def load_client_data(partition_id: int):
    
    # Load data for the given partition
    config_file_path = 'configs/center_' + str(partition_id) + '/config.yaml'
    loader = DataLoader(config_file_path)
    
    return loader
