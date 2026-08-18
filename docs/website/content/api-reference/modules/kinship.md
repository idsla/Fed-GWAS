---
title: KING / Kinship
---

# KING / Kinship

This module implements **KING-based relatedness screening**. Flower stage `iterative_king` drives kinship estimation: clients send sample-based PLINK 1.9 chunks; the server reconstructs, merges, and runs PLINK 2 `--make-king-table`.

## Client flow

`handle_iterative_king` checks whether the client already has `chunk_data`. If not, it recreates sample chunks:

```python
partition_config = {
    "chunk_size": chunk_size,
    "partition_by": "samples",
}
client.partition_data(partition_config)
```

Each call sends one chunk and increments `current_chunk_idx`. After all chunks are sent, the client finalizes local KING accumulators and can filter samples above `king_threshold`.

## Server flow

`run_server_king`:

1. Decodes each received chunk array.
2. Rebuilds temporary `.bed/.bim/.fam` files from size metadata.
3. Merges PLINK datasets.
4. Runs PLINK `--het`.
5. Runs PLINK2 `--make-king-table`.
6. Parses `.kin0` output and returns sample pairs, partial phi values, and heterozygosity-derived information.

## Outputs

Server-side temporary files and KING outputs are written into a session directory under the strategy `output_dir`. Clients also persist `king_accumulator_state.pkl` and `king_results_<client_id>.txt` under their configured log directories.
