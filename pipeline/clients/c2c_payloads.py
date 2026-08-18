from __future__ import annotations

import numpy as np


def pack_typed_payload_uint8(kind: str, arr: np.ndarray) -> np.ndarray:
    """
    Pack a typed payload into a uint8 ndarray, suitable for encrypting with ClientToClientMessenger.

    Format:
      b"FGQP" + kind_len(u8) + kind_ascii + encoded_ndarray_bytes
    """
    from pipeline.clients.client_to_client import _encode_ndarray  # internal helper

    kind_b = kind.encode("ascii")
    if len(kind_b) > 255:
        raise ValueError("kind too long")

    body = _encode_ndarray(arr)
    blob = b"FGQP" + bytes([len(kind_b)]) + kind_b + body
    return np.frombuffer(blob, dtype=np.uint8)


def unpack_typed_payload_uint8(arr_u8: np.ndarray) -> tuple[str, np.ndarray]:
    """Inverse of pack_typed_payload_uint8()."""
    from pipeline.clients.client_to_client import _decode_ndarray  # internal helper

    blob = arr_u8.tobytes()
    if len(blob) < 5 or blob[:4] != b"FGQP":
        raise ValueError("bad payload magic")

    kind_len = blob[4]
    if len(blob) < 5 + kind_len:
        raise ValueError("bad payload kind length")

    kind = blob[5 : 5 + kind_len].decode("ascii")
    body = blob[5 + kind_len :]
    return kind, _decode_ndarray(body)
