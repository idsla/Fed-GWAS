# server/prg_masking.py

import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from Crypto.Protocol.KDF import HKDF
from Crypto.PublicKey import ECC
from Crypto.Hash import SHA256
from Crypto.Protocol.DH import key_agreement

class PRGMaskingAggregator:
    """
    PRG-MASKING secure aggregation based on CCS 2017 paper.
    Implements pairwise masking with Diffie-Hellman key exchange.
    """
    
    def __init__(self, num_clients: int):
        self.num_clients = num_clients
        # Use P-256 curve for secure key exchange
        self.curve = 'p256'
        
        # Store client public keys (as ECC key objects)
        self.client_public_keys: Dict[int, ECC.EccKey] = {}
        
        # State management
        self.key_exchange_complete = False
        
    def get_curve_params(self) -> Dict[str, str]:
        """Get curve parameters for key exchange."""
        return {"curve": self.curve}
    
    def add_client_public_key(self, client_id: int, public_key_pem: str):
        """Store client's public key."""
        try:
            # Import PEM format public key
            public_key = ECC.import_key(public_key_pem)
            self.client_public_keys[client_id] = public_key
            
            # Check if key exchange is complete
            if len(self.client_public_keys) == self.num_clients:
                self.key_exchange_complete = True
                print(f"[PRG-MASKING] Key exchange complete with {self.num_clients} clients")
        except Exception as e:
            print(f"[PRG-MASKING] Error importing public key for client {client_id}: {e}")
    
    def is_key_exchange_complete(self) -> bool:
        """Check if all clients have submitted public keys."""
        return self.key_exchange_complete
    
    def get_all_public_keys(self) -> Dict[str, str]:
        """Get all client public keys for distribution."""
        return {str(k): v.export_key(format='PEM') for k, v in self.client_public_keys.items()}
    
    def generate_shared_secrets(self, client_id: int, private_key: ECC.EccKey) -> Dict[int, bytes]:
        """
        Generate shared secrets with all other clients.
        This is used by clients to compute pairwise masks.
        """
        shared_secrets = {}
        salt = b'fedgwas-prg-masking-salt'
        
        def kdf(shared_material: bytes) -> bytes:
            return HKDF(
                master=shared_material,
                key_len=32,
                salt=salt,
                hashmod=SHA256
            )
        
        for peer_id in range(self.num_clients):
            if peer_id == client_id:
                continue
                
            peer_pub_key = self.client_public_keys[peer_id]
            
            # Generate shared key material using ECC key agreement
            shared_key_material = key_agreement(
                static_priv=private_key,
                static_pub=peer_pub_key,
                kdf=kdf
            )
            
            shared_secrets[peer_id] = shared_key_material
            
        return shared_secrets

class ClientMaskingHelper:
    """
    Helper class for client-side PRG masking operations.
    """
    
    def __init__(self, client_id: int, num_clients: int):
        self.client_id = client_id
        self.num_clients = num_clients
        self.private_key = None
        self.public_key_pem = None
        self.shared_secrets: Dict[int, bytes] = {}
    
    def generate_ecc_keypair(self, curve_params: Dict[str, str]) -> str:
        """Generate ECC key pair and return public key PEM string."""
        curve = curve_params.get("curve", "p256")
        
        # Generate ECC key pair
        self.private_key = ECC.generate(curve=curve)
        
        # Export public key in PEM format
        self.public_key_pem = self.private_key.public_key().export_key(format='PEM')
        
        return self.public_key_pem
    
    def compute_shared_secrets(self, all_public_keys: Dict[str, str], curve_params: Dict[str, str]):
        """Compute shared secrets with all other clients."""
        salt = b'fedgwas-prg-masking-salt'
        
        def kdf(shared_material: bytes) -> bytes:
            return HKDF(
                master=shared_material,
                key_len=32,
                salt=salt,
                hashmod=SHA256
            )
        
        for peer_id_str, peer_pub_key_pem in all_public_keys.items():
            peer_id = int(peer_id_str)
            if peer_id == self.client_id:
                continue
            
            try:
                # Import peer's public key from PEM
                peer_pub_key = ECC.import_key(peer_pub_key_pem)
                
                # Generate shared key material using ECC key agreement
                shared_key_material = key_agreement(
                    static_priv=self.private_key,
                    static_pub=peer_pub_key,
                    kdf=kdf
                )
                
                self.shared_secrets[peer_id] = shared_key_material
            except Exception as e:
                print(f"[Client {self.client_id}] Error computing shared secret with peer {peer_id}: {e}")
    
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