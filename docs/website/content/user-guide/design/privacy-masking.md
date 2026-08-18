---
title: Privacy and Masking
slug: /design/privacy-masking
---

# Privacy and Masking

FedGWAS uses encryption, shuffling, anonymization, and lightweight secret-sharing to protect sensitive information during federated execution. The coordinating server relays selected client-to-client payloads and does not decrypt them. These mechanisms provide practical protection against direct information leakage; they are not a complete privacy proof for every deployment threat model.

The current code implements those paper-level mechanisms as ECC-backed encrypted envelopes, privacy-preserving tokens for association filtering, and deterministic identifier anonymization for chunked PLINK data. Helper modules such as `prg_masking.py` support key exchange for the encrypted relay; they are not a separate “secure aggregation” protocol.

## Encrypted relay

`pipeline/server/prg_masking.py` tracks public keys and exposes curve parameters. `pipeline/clients/client_to_client.py` uses client private keys and peer public keys to build encrypted envelopes that the server forwards but cannot decrypt.

The flow is:

1. During `key_exchange` (initialization), each client generates an ECC keypair and sends a PEM-encoded public key.
2. The server stores public keys and distributes the full public-key map during later stages.
3. Each client computes shared secrets with peers (lightweight secret-sharing of seed material).
4. Clients encrypt seed messages, QC payloads, and relatedness-map payloads for peer recipients.
5. The server only relays ciphertext arrays through Flower `Parameters`.
6. Recipient clients decrypt the envelopes and perform the required local aggregation or mapping.

The same relay pattern is used for seed synchronization and federated quality control. Local association filtering sends privacy-preserving (salted) tokens instead of raw SNP IDs; clients compute the intersection locally during `local_lr_filter_response`.

## ID anonymization

`BaseGWASClient` contains helpers for anonymizing samples and SNPs before chunk transport. Combined with shuffling of chunk order under the shared global seed, this reduces linkage of identifiers across chunks or federated rounds.

| Data | Method |
| --- | --- |
| Sample FID/IID | Add `sample_offset` to numeric IDs; use a hash fallback for non-numeric IDs. |
| SNP ID | Compute a deterministic SHA-256 based ID from `global_seed` and the original SNP ID. |

The client persists chunk map files under its intermediate directory so returned anonymous IDs can be mapped back to local identifiers when possible, including after Flower actor restarts.

## Operational notes

- `num_clients` must match the actual number of clients so each client encrypts messages for the right recipients.
- Simulation mode may recreate Flower actors across rounds; the client stores configuration in Flower state to mitigate that.
- Encrypted relay, shuffling, anonymization, and tokenization reduce server visibility. Production use should include a separate security review.
