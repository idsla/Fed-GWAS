---
title: Troubleshooting
slug: /troubleshooting
---

# Troubleshooting

## PLINK is not found

Typical symptom:

```text
PLINK command failed
```

Check the following:

- `plink` or `plink.exe` is available on `PATH`.
- On Windows, the project `bin/` directory contains `plink.exe` or is added to `PATH`.
- `input_data.path` points to a PLINK prefix and does not include `.bed`.

## Missing curve parameters during `sync`

Typical symptom:

```text
Missing curve params in sync stage
```

Check the following:

- The server completed the `key_exchange` stage and collected all client public keys.
- The number of Flower clients matches the masking helper's expected client count.
- In simulation mode, inspect server warnings for relay stages that received no client results; `strategy_strict.py` retries these stages instead of advancing on failures.

## Client logs are empty

Client logs are written under each center's configured `output.log_dir` with per-client filenames.

Check the following:

- `output.log_dir` exists or can be created.
- Client code logs through `self.logger`.
- The client did not fail before `LoggerManager.get_logger(...)` was called.

## Deployment port errors

Expected local-deployment ports:

| Port | Role |
| --- | --- |
| `9092` | Fleet API for SuperNodes |
| `9093` | Exec API for `flwr run` |
| `9094+` | ClientAppIo API for SuperNodes |

Connection refused errors usually mean SuperLink or the relevant SuperNode is not running.
