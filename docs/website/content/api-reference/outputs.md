---
title: Output Files
---

# Output Files

FedGWAS writes outputs to the paths declared in each experiment's `configs/server/config.yaml` and `configs/center_x/config.yaml`.

## Result tree

The shipped tiny correctness config writes to `experiments/correctness/tiny_even/results_2/`:

```text
results_2/
  server/
    logs/
      server_log.txt
      stage_metrics.csv
    intermediate/
      king_session_*/
      lr_session_*/
  center_1/
    logs/
      client_0_log.txt
      stage_metrics.csv
      client_metrics.csv
    intermediate/
      run_<run_id>/
  center_2/
    logs/
    intermediate/
  stage_metrics.csv
  client_metrics.csv
  retention_manifest.json
```

Exact filenames can vary by stage and retention tier.

## Important files

| File | Producer | Purpose |
| --- | --- | --- |
| `server/logs/server_log.txt` | Server strategy | Stage progress, relay behavior, KING/LR execution, completion status. |
| `center_x/logs/client_*_log.txt` | Client app | Per-client stage progress and PLINK command details. |
| `stage_metrics.csv` | Monitoring runtime | Stage duration and status metrics. |
| `client_metrics.csv` | Client monitoring runtime | Client-side resource and network samples when enabled. |
| `king_results_<client_id>.txt` | Client KING stage | Local accumulated kinship output after de-anonymization where possible. |
| `lr_results_latest.assoc.logistic` | Server LR aggregator | Latest server-side PLINK logistic output copy. |
| `retention_manifest.json` | Retention utility | Files kept or pruned at completion. |
| `qc_report.md`, `lr_report.md`, `evaluation_report.md` | Evaluation tools | Post-run agreement reports against centralized baseline. |

## Retention

Experiment `config.yaml` can set:

```yaml
retention:
  tier: standard
  auto_apply_on_complete: true
```

Retention is applied when the server reaches `done`. Use `experiments/tools/apply_run_retention.py --dry-run` to preview manual pruning.
