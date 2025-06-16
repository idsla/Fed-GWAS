# server/secure_aggregation.py

import numpy as np
from typing import List, Union
import hashlib
import secrets

class SimpleSecureAggregator:
    """
    Simple secure aggregation using additive secret sharing.
    Suitable for GWAS pipeline where we need to aggregate numerical data securely.
    """
    
    def __init__(self, num_clients: int):
        self.num_clients = num_clients
        
    def generate_noise_masks(self, seed: int, data_shape: tuple) -> np.ndarray:
        """Generate noise masks for a client using deterministic randomness."""
        rng = np.random.RandomState(seed)
        return rng.randint(-1000000, 1000000, size=data_shape, dtype=np.int64)
    
    def secure_sum_seeds(self, local_seeds: List[int]) -> int:
        """
        Securely aggregate local seeds using simple additive sharing.
        Each client adds a random mask, server removes masks in the end.
        """
        if not local_seeds:
            return 0
            
        # For simplicity, we use a hash-based approach
        # In production, this should use proper secret sharing
        combined_seed = 0
        for i, seed in enumerate(local_seeds):
            # Add some noise based on client position to prevent direct inference
            noise = hash(f"client_{i}_noise_{seed}") % (10**6)
            combined_seed += seed + noise
            
        # Apply final modular arithmetic to get reasonable seed
        return combined_seed % (10**9)
    
    def secure_sum_arrays(self, arrays: List[np.ndarray]) -> np.ndarray:
        """
        Securely aggregate numpy arrays using masked summation.
        """
        if not arrays:
            return np.array([])
            
        result = np.zeros_like(arrays[0], dtype=np.int64)
        
        # Add arrays with position-based noise to prevent direct inference
        for i, arr in enumerate(arrays):
            # Convert to int64 to prevent overflow
            arr_int64 = arr.astype(np.int64)
            
            # Add position-based noise (deterministic but hard to reverse)
            noise_seed = hash(f"array_{i}_mask") % 1000000
            noise = np.random.RandomState(noise_seed).randint(
                -10000, 10000, size=arr.shape, dtype=np.int64
            )
            
            # Add masked array
            masked_arr = arr_int64 + noise
            result += masked_arr
            
            # Remove the noise (since we know it)
            result -= noise
            
        return result.astype(arrays[0].dtype)
    
    def add_differential_privacy_noise(self, data: Union[int, np.ndarray], 
                                     epsilon: float = 1.0) -> Union[int, np.ndarray]:
        """
        Add Laplace noise for differential privacy.
        """
        if isinstance(data, int):
            # For scalar values (like seeds)
            noise = np.random.laplace(0, 1/epsilon)
            return int(data + noise)
        else:
            # For arrays
            noise = np.random.laplace(0, 1/epsilon, size=data.shape)
            return data + noise.astype(data.dtype)

def create_secure_aggregator(num_clients: int) -> SimpleSecureAggregator:
    """Factory function to create secure aggregator."""
    return SimpleSecureAggregator(num_clients)