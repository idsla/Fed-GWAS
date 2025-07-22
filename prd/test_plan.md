# Testing Plan for Federated GWAS Pipeline

This document outlines how to validate and benchmark the federated GWAS pipeline under different deployment and data‐size scenarios, **given that late‐join clients are not supported**. Any client must start at the very first stage (`local_qc`) or it will be excluded from all downstream stages.

---

## 1. Objectives

1. **Correctness**: Verify that each stage (Local QC, Global QC, Sync, Iterative KING, Local LR, Iterative LR) produces the expected outputs.  
2. **Robustness**: Ensure the pipeline tolerates clients joining only at the beginning and gracefully ignores any attempt to rejoin mid-pipeline.  
3. **Performance**: Measure wall-clock time, memory usage, and network I/O under varying numbers of clients, machines, and data sizes.  

---

## 2. Test Scenarios

| Scenario | Machines | Clients per Machine | Total Clients | Notes                                             |
|:--------:|:--------:|:-------------------:|:-------------:|:-------------------------------------------------:|
| **A**    | 1        | 1                   | 1             | Baseline, single-client end-to-end                |
| **B**    | 1        | 3                   | 3             | All three clients start together, shared host     |
| **C**    | 3        | 1                   | 3             | One client per machine, concurrent participation |

*Late-join is disallowed: any client not present at QC (or that opts out early) will be filtered out and not participate in later stages.*

---

## 3. Dataset Sizes & Partitioning

| Size    | #Samples | #SNPs    | Approx. BED Size | Use Case               |
|:-------:|:--------:|:--------:|:----------------:|:----------------------:|
| **Small**  | 100      | 10 000   | ~50 MB           | Quick functional tests |
| **Medium** | 1 000    | 100 000  | ~500 MB          | Realistic development  |
| **Large**  | 10 000   | 1 000 000| ~5 GB            | Scalability benchmarks |

Test each scenario (A/B/C) at **Small**, **Medium**, and **Large** scales.

**Another approach is to investigate the exisiting framework and machtch the scale in their implementation**

---

## 4. Participation Patterns

Since in current version, **late joins are not supported**, clients must either:

- **Fully participate** from `local_qc` through `iterative_lr`, or  
- **Opt out entirely** at certain stage (and then they exit via `sys.exit(0)`).

Any client that opts out or fails to participate in the very first sync will be excluded from all downstream rounds.

---

## 5. Test Matrix

| Scenario | Data Size | Participation Pattern             | Expected Outcome                                          |
|:--------:|:---------:|:---------------------------------:|:---------------------------------------------------------:|
| A1       | Small     | Full participation (1 client)     | End-to-end success                                        |
| A2       | Medium    | Full participation                | Performance scaling                                       |
| A3       | Large     | Full participation                | Resource usage profiling                                  |
| B1       | Small     | 2 clients drop before QC          | Remaining 1 client completes; others exit                 |
| B2       | Medium    | 1 client drops after QC           | Two complete end-to-end; dropped client never re-joins     |
| B3       | Large     | All participate initially         | All three complete; measure aggregate performance         |
| C1       | Small     | 1 client drops during KING        | Two continue KING and LR; dropped client excluded downstream |
| C2       | Medium    | 1 client drops during LR          | Two complete LR; dropped client excluded                  |
| C3       | Large     | All complete                      | Full three-client run; benchmark at scale                 |

---

## 6. Test Steps

For each test case:

1. **Environment Setup**  
   - Install PLINK, Python dependencies, Flower.  
   - Populate `client/config.yaml` with appropriate `participation` flags and directories.  
   - Place input data under `input_data.path`.

2. **Launch Server**  
   ```bash
   cd server/
   python main_server.py
   
3. **Steps**
   - Clients write to logs/iteration_log.txt (timestamps, stage, actions, exits).
   - Server logs participants_per_stage and failures_per_stage.

4. **Collect Metrics**
   - Timing: record start/end of each stage per client.
   - Memory:
   - Network I/O: 

5. **Validate Outputs**
   - QC: Check that excluded SNP lists match manual runs.
   - KING: Verify kinship coefficients sum correctly across chunks.
   - LR: Compare top p-values list against a centralized run.

6. **Document Results**
   - For each run, record in a spreadsheet:
   - Scenario ID, Data Size, Dynamic Pattern
   - Wall-clock time per stage (avg client & server)
   - Memory peak
   - Total bytes sent/received
   - Any failures/exceptions

7. **Reporting**
   - Create a summary report with tables and graphs:
   - Stage runtimes vs. data size
   - Client participation rates vs. pipeline progression
   - Error/failure counts by stage and cause