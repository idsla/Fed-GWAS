from __future__ import annotations

import json
from typing import Dict


def get_all_public_keys(config: dict) -> Dict[str, str]:
    """
    Read all_public_keys from Flower config without relying on dict-typed Scalars.

    Strategy may pass either:
    - `all_public_keys` as a dict (rare, depends on Flower Scalar support)
    - `all_public_keys_json` as a JSON string (preferred)
    """
    apk = config.get("all_public_keys")
    if isinstance(apk, dict):
        return apk

    apk_json = config.get("all_public_keys_json")
    if isinstance(apk_json, str) and apk_json:
        try:
            parsed = json.loads(apk_json)
            return parsed if isinstance(parsed, dict) else {}
        except Exception:
            return {}

    return {}


