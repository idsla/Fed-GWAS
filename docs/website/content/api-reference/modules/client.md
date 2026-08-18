---
title: Client Modules
---

# Client Modules

The client entry point is `pipeline/client_app.py`. Flower calls `client_fn(context)`, which creates a `FedLRClient` for the current partition.

## Initialization path

`FedLRClient.__init__` performs the following work:

1. Reads the partition id from the Flower `Context`.
2. Resolves the center configuration path. In simulation mode, this is derived from `config_path` and `partition-id`.
3. Loads the YAML file with `DataLoader`.
4. Converts or binds the input data to a PLINK prefix.
5. Stores configuration records in Flower client state on first initialization.
6. Restores configuration records from state on later rounds.
7. Initializes `BaseGWASClient`, thresholds, participation flags, and the per-client logger.

## Stage routing

`fit(parameters, config)` reads `config["stage"]` and delegates to the matching branch:

| Stage | Implementation |
| --- | --- |
| `key_exchange` | `pipeline/client_app.py` with helpers from `prg_masking.py` |
| `sync` / `sync_response` | `pipeline/clients/seed_sync.py` and `pipeline/clients/client_to_client.py` |
| `local_qc` | `pipeline/clients/local_qc.py` |
| `global_qc` / `global_qc_response` | `pipeline/clients/local_qc.py`, `pipeline/clients/c2c_payloads.py`, and `pipeline/clients/client_qc_aggregator.py` |
| `local_lr` / `local_lr_filter_response` | `pipeline/clients/local_qc.py` and `pipeline/clients/lr_privacy.py` |
| `init_chunks` / `init_chunks_lr` | `BaseGWASClient.partition_data` |
| `iterative_king` | `pipeline/clients/iterative_king.py` |
| `iterative_lr` / `iterative_lr_response` | `pipeline/clients/iterative_lr.py` |

## Chunking behavior

`BaseGWASClient.partition_data` supports three partition modes:

- `samples`: reads `.fam`, shuffles sample IDs with `global_seed`, writes a PLINK `--keep` file, and creates sample chunks.
- `snps`: reads `.bim`, shuffles SNP IDs with `global_seed`, writes a PLINK `--extract` file, and creates SNP chunks.
- `both`: currently implemented as a simplified single-chunk path.

Each chunk is rebuilt as `.bed/.bim/.fam`, encoded into a numpy array, sent to the server, and then removed from the local intermediate directory. Intermediate directories can be scoped to `run_<run_id>` when `run_id` or `run_tag` is configured.
