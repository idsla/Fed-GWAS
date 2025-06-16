# server/prg_masking.py

import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from Cryptodome.Protocol.KDF import HKDF
from Cryptodome.PublicKey import DH
from Cryptodome.Hash import SHA256

class PRGMaskingAggregator:
    """
    PRG-MASKING secure aggregation based on CCS 2017 paper.
    Implements pairwise masking with Diffie-Hellman key exchange.
    """
    
    def __init__(self, num_clients: int):
        self.num_clients = num_clients
        # Generate global DH parameters (use 1024 for demo, 2048+ for production)
        self.dh_obj = DH.generate(1024)
        self.dh_params = {
            "dh_p": str(self.dh_obj.p), 
            "dh_g": str(self.dh_obj.g)
        }
        
        # Store client public keys
        self.client_public_keys: Dict[int, str] = {}
        
        # State management
        self.key_exchange_complete = False
        
    def get_dh_params(self) -> Dict[str, str]:
        """Get DH parameters for key exchange."""
        return self.dh_params.copy()
    
    def add_client_public_key(self, client_id: int, public_key_str: str):
        """Store client's public key."""
        self.client_public_keys[client_id] = public_key_str
        
        # Check if key exchange is complete
        if len(self.client_public_keys) == self.num_clients:
            self.key_exchange_complete = True
            print(f"[PRG-MASKING] Key exchange complete with {self.num_clients} clients")
    
    def is_key_exchange_complete(self) -> bool:
        """Check if all clients have submitted public keys."""
        return self.key_exchange_complete
    
    def get_all_public_keys(self) -> Dict[str, str]:
        """Get all client public keys for distribution."""
        return {str(k): v for k, v in self.client_public_keys.items()}
    
    def generate_shared_secrets(self, client_id: int, private_key: int) -> Dict[int, bytes]:
        """
        Generate shared secrets with all other clients.
        This is used by clients to compute pairwise masks.
        """
        shared_secrets = {}
        salt = b'fedgwas-prg-masking-salt'
        
        for peer_id in range(self.num_clients):
            if peer_id == client_id:
                continue
                
            peer_pub_key_str = self.client_public_keys[peer_id]
            peer_pub_key = DH.construct((int(peer_pub_key_str), self.dh_obj.p, self.dh_obj.g))
            
            # Generate shared key material using DH
            shared_key_material = self.dh_obj.exchange(private_key, peer_pub_key)
            
            # Derive final shared key using HKDF
            derived_key = HKDF(
                master=shared_key_material, 
                key_len=32, 
                salt=salt, 
                hashmod=SHA256
            )
            shared_secrets[peer_id] = derived_key
            
        return shared_secrets

class ClientMaskingHelper:
    """
    Helper class for client-side PRG masking operations.
    """
    
    def __init__(self, client_id: int, num_clients: int):
        self.client_id = client_id
        self.num_clients = num_clients
        self.private_key = None
        self.public_key_str = None
        self.shared_secrets: Dict[int, bytes] = {}
    
    def generate_dh_keypair(self, dh_params: Dict[str, str]) -> str:
        """Generate DH key pair and return public key string."""
        # Reconstruct DH object from parameters
        dh_obj = DH.construct((int(dh_params["dh_p"]), int(dh_params["dh_g"])))
        
        # Generate private key
        self.private_key = dh_obj.generate_private_key()
        
        # Generate public key
        public_key = dh_obj.generate_public_key(self.private_key)
        self.public_key_str = str(public_key)
        
        return self.public_key_str
    
    def compute_shared_secrets(self, all_public_keys: Dict[str, str], dh_params: Dict[str, str]):
        """Compute shared secrets with all other clients."""
        dh_obj = DH.construct((int(dh_params["dh_p"]), int(dh_params["dh_g"])))
        salt = b'fedgwas-prg-masking-salt'
        
        for peer_id_str, peer_pub_key_str in all_public_keys.items():
            peer_id = int(peer_id_str)
            if peer_id == self.client_id:
                continue
            
            peer_pub_key = DH.construct((int(peer_pub_key_str), dh_obj.p, dh_obj.g))
            
            # Generate shared key material
            shared_key_material = dh_obj.exchange(self.private_key, peer_pub_key)
            
            # Derive shared secret
            derived_key = HKDF(
                master=shared_key_material,
                key_len=32,
                salt=salt,
                hashmod=SHA256
            )
            self.shared_secrets[peer_id] = derived_key
    
    def mask_data(self, data: np.ndarray) -> np.ndarray:
        """
        Apply PRG masking to data using pairwise shared secrets.
        
        Formula: y_i = x_i + sum_{j>i} p_ij - sum_{j<i} p_ji
        where p_ij is generated from shared secret s_ij
        """
        if isinstance(data, (int, float)):
            # Handle scalar data (like seeds)
            data = np.array([data], dtype=np.float64)
            is_scalar = True
        else:
            # Handle array data (like QC arrays)
            data = data.astype(np.float64)
            is_scalar = False
        
        masking_vector = np.zeros_like(data, dtype=np.float64)
        
        for peer_id in range(self.num_clients):
            if peer_id == self.client_id:
                continue
            
            shared_secret = self.shared_secrets[peer_id]
            
            # Generate pseudorandom vector using shared secret as seed
            prng = np.random.RandomState(seed=np.frombuffer(shared_secret[:4], dtype=np.uint32)[0])
            p_vector = prng.normal(0, 1, size=data.shape)  # Use normal distribution
            
            # Apply masking based on client ID ordering
            if self.client_id < peer_id:
                masking_vector += p_vector  # Add for j > i
            else:
                masking_vector -= p_vector  # Subtract for j < i
        
        masked_data = data + masking_vector
        
        if is_scalar:
            return masked_data[0]
        else:
            return masked_data

def create_prg_masking_aggregator(num_clients: int) -> PRGMaskingAggregator:
    """Factory function to create PRG masking aggregator."""
    return PRGMaskingAggregator(num_clients)

def create_client_masking_helper(client_id: int, num_clients: int) -> ClientMaskingHelper:
    """Factory function to create client masking helper."""
    return ClientMaskingHelper(client_id, num_clients)