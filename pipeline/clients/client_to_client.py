"""
Client-to-Client Communication via Encrypted Messages

This module implements encrypted client-to-client communication where:
1. Clients encrypt messages for specific recipients using the recipient's public key
2. Server forwards encrypted messages but cannot decrypt them (doesn't have private keys)
3. Recipients decrypt messages using their private keys

This enables true client-to-client communication even in a server-mediated architecture.
"""

from __future__ import annotations

import json
import struct
from typing import Dict

import numpy as np
from Crypto.Cipher import AES
from Crypto.Hash import SHA256
from Crypto.Protocol.DH import key_agreement
from Crypto.Protocol.KDF import HKDF
from Crypto.PublicKey import ECC
from Crypto.Random import get_random_bytes
from Crypto.Util.Padding import pad, unpad


def _encode_ndarray(arr: np.ndarray) -> bytes:
    """Serialize a numeric ndarray without pickle."""
    if not isinstance(arr, np.ndarray):
        arr = np.asarray(arr)
    header = {"dtype": str(arr.dtype), "shape": list(arr.shape), "order": "C"}
    header_bytes = json.dumps(header, separators=(",", ":")).encode("utf-8")
    payload = arr.tobytes(order="C")
    return struct.pack("!I", len(header_bytes)) + header_bytes + payload


def _decode_ndarray(blob: bytes) -> np.ndarray:
    """Deserialize ndarray from _encode_ndarray."""
    if len(blob) < 4:
        raise ValueError("Invalid ndarray blob (too short)")
    (header_len,) = struct.unpack("!I", blob[:4])
    header_start = 4
    header_end = header_start + header_len
    header = json.loads(blob[header_start:header_end].decode("utf-8"))
    dtype = np.dtype(header["dtype"])
    shape = tuple(header["shape"])
    payload = blob[header_end:]
    arr = np.frombuffer(payload, dtype=dtype)
    return arr.reshape(shape)


def pack_envelope(sender_id: int, recipient_id: int, iv: bytes, encrypted_data: bytes) -> bytes:
    """
    Binary envelope format (no pickle):
      magic(4) +
      sender_len(uint16) + sender_id_ascii +
      recip_len(uint16) + recipient_id_ascii +
      iv(16) +
      cipher_len(uint32) + ciphertext

    IDs can exceed int64; encode as ASCII to avoid overflow.
    """
    if len(iv) != 16:
        raise ValueError("IV must be 16 bytes")
    magic = b"FGCC"
    sender_b = str(int(sender_id)).encode("ascii")
    recip_b = str(int(recipient_id)).encode("ascii")
    if len(sender_b) > 65535 or len(recip_b) > 65535:
        raise ValueError("sender/recipient id too long to encode")
    header = struct.pack("!4sH", magic, len(sender_b)) + sender_b
    header += struct.pack("!H", len(recip_b)) + recip_b
    header += iv
    header += struct.pack("!I", len(encrypted_data))
    return header + encrypted_data


def unpack_envelope(envelope: bytes) -> Dict[str, object]:
    """Inverse of pack_envelope."""
    if len(envelope) < 4 + 2 + 2 + 16 + 4:
        raise ValueError("Invalid envelope (too short)")
    off = 0
    magic, sender_len = struct.unpack_from("!4sH", envelope, off)
    off += 6
    if magic != b"FGCC":
        raise ValueError("Invalid envelope (bad magic)")
    sender_b = envelope[off : off + sender_len]
    off += sender_len
    (recip_len,) = struct.unpack_from("!H", envelope, off)
    off += 2
    recip_b = envelope[off : off + recip_len]
    off += recip_len
    iv = envelope[off : off + 16]
    off += 16
    (clen,) = struct.unpack_from("!I", envelope, off)
    off += 4
    ciphertext = envelope[off : off + clen]
    if len(ciphertext) != clen:
        raise ValueError("Invalid envelope (ciphertext length mismatch)")
    return {
        "sender_id": int(sender_b.decode("ascii")),
        "recipient_id": int(recip_b.decode("ascii")),
        "iv": iv,
        "encrypted_data": ciphertext,
    }


class ClientToClientMessenger:
    """Encrypt/decrypt messages routed via server without server decryption ability."""

    def __init__(self, client_id: int, num_clients: int, private_key: ECC.EccKey):
        self.client_id = client_id
        self.num_clients = num_clients
        self.private_key = private_key
        self.shared_keys: Dict[int, bytes] = {}

    def _derive_aes_key(self, peer_id: int, peer_public_key: ECC.EccKey) -> bytes:
        if peer_id in self.shared_keys:
            return self.shared_keys[peer_id]

        salt = b"fedgwas-client-to-client-salt"

        def kdf(shared_material: bytes) -> bytes:
            return HKDF(master=shared_material, key_len=32, salt=salt, hashmod=SHA256)

        aes_key = key_agreement(static_priv=self.private_key, static_pub=peer_public_key, kdf=kdf)
        self.shared_keys[peer_id] = aes_key
        return aes_key

    def encrypt_for_recipient(self, data: np.ndarray, recipient_id: int, recipient_public_key_pem: str) -> bytes:
        recipient_public_key = ECC.import_key(recipient_public_key_pem)
        aes_key = self._derive_aes_key(recipient_id, recipient_public_key)
        data_bytes = _encode_ndarray(data)
        iv = get_random_bytes(16)
        cipher = AES.new(aes_key, AES.MODE_CBC, iv)
        encrypted = cipher.encrypt(pad(data_bytes, AES.block_size))
        return pack_envelope(self.client_id, recipient_id, iv, encrypted)

    def decrypt_from_sender(self, encrypted_message: bytes, sender_id: int, sender_public_key_pem: str) -> np.ndarray:
        metadata = unpack_envelope(encrypted_message)
        sender_public_key = ECC.import_key(sender_public_key_pem)
        aes_key = self._derive_aes_key(sender_id, sender_public_key)
        iv = metadata["iv"]
        encrypted_data = metadata["encrypted_data"]
        cipher = AES.new(aes_key, AES.MODE_CBC, iv)
        decrypted = unpad(cipher.decrypt(encrypted_data), AES.block_size)
        return _decode_ndarray(decrypted)


def create_client_messenger(client_id: int, num_clients: int, private_key: ECC.EccKey) -> ClientToClientMessenger:
    return ClientToClientMessenger(client_id, num_clients, private_key)

