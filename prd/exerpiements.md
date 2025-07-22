## 1. Simulated Datasets

We generate synthetic PLINK‐format datasets at four scales. Each scale is partitioned horizontally (by samples) across a fixed number of centers.

### 1.1. Data Generation by Scale

1. **Tiny (Sanity Check)**
   - **Samples:** 500  
   - **SNPs:** 5 000  
   - **Centers (Clients):** 2  
   - **Purpose:** Quick correctness checks and dev‐time regression.

2. **Small (Functional Test)**
   - **Samples:** 5 000  
   - **SNPs:** 50 000  
   - **Centers (Clients):** 3  
   - **Purpose:** End-to-end functionality, basic performance profiling (minutes to hours).

3. **Medium (Realistic Development)**
   - **Samples:** 50 000  
   - **SNPs:** 500 000  
   - **Centers (Clients):** 5  
   - **Purpose:** Early performance tuning (runtime ~hours), memory profiling.

4. **Large (UKB-Scale Proxy)**
   - **Samples:** 275 812 (≈ UK Biobank)  
   - **SNPs:** 98 000 000 (≈ imputed genotypes)  
   - **Centers (Clients):** 7  
   - **Purpose:** Scalability benchmark, compare to SF‐GWAS published runtimes (QC ~4.5 h, PCA ~44 h, LR ~77.8 h, total ~5.3 days).  
   - **Goal:** Demonstrate that our federated pipeline can approach or improve on those runtime figures; capture per‐stage wall‐clock and resource usage.


**Note:** Validate that for the “Large” scale, the total combined PLINK binary files approach ~5 TB (estimated from UKB Imputed). If disk constraints exist, generate a subset of ~10 M–20 M SNPs to mimic large‐scale behavior (e.g. 1/5 scale) and extrapolate.

#### 1.1.1. Allele Frequency & Genotype Simulation

- Draw each SNP’s minor‐allele frequency (MAF) from Uniform(0.01, 0.5).
- For each sample–SNP pair, sample genotype ∈ {0, 1, 2} via Binomial(2, p_snp).
- Introduce ~1–3 % missing calls at random to simulate real‐world missingness.

#### 1.1.2. Phenotype Simulation

- **Continuous trait:**  
  - Randomly select 1 % of SNPs as “causal,” assign effect sizes β ∼ N(0, 0.1).  
  - Phenotype Y = Σ(βᵢ·Gᵢ) + ε, ε ∼ N(0, 1).  

- **Binary trait (optional):**  
  - Use a logistic model on the same set of causal SNPs to generate case/control labels.

#### 1.1.3. Horizontal Partitioning by Samples

1. **Even Sample Split**  
   - Divide total samples equally across centers.  
   - Each center has the same phenotype distribution (balanced cases/controls).

2. **Uneven Sample Sizes**  
   - Assign different sample counts per center.  
   - Phenotype distribution remains balanced (50/50) in each center.

3. **Uneven Size + Phenotype Distribution**  
  e.g.,
   - Center A: large sample count (e.g., 100 k) with balanced phenotype (≈50 % cases, 50 % controls).  
   - Center B: medium sample count (e.g., 50 k) with moderately skewed phenotype (e.g., 30 % cases, 70 % controls).  
   - Center C: smaller sample count (e.g., 25 k) with highly skewed phenotype (e.g., 10 % cases, 90 % controls).  
   - Remaining centers: similar logic if more than three.  
   - All centers hold the same full SNP set; only sample counts and phenotype ratios differ.

**Implementation notes:**  
1. Generate a single full‐scale dataset codebook of sample IDs (with assigned phenotypes).  
2. Partition the sample IDs into center‐specific lists according to the chosen strategy (e.g., uneven sizes + phenotype ratios).  
3. For each center, run PLINK with `--keep` on that center’s sample list to produce a local `.bed/.bim/.fam`.  
4. All centers share an identical `.bim` and `.fam` schema (same SNP IDs, order, map).

<!-- 

- **Implementation:**  
  1. Generate full PLINK dataset (one `.bed/.bim/.fam`) at the desired scale.  
  2. Compute the list of sample IDs.  
  3. For each center, write a sample‐keep file listing that center’s sample IDs, then run:  
     ```bash
     plink --bfile full_dataset_prefix \
           --keep centerX_sample_list.txt \
           --make-bed \
           --out centerX_dataset_prefix
     ```  
  4. Repeat for each center. Now each centerX has `.bed/.bim/.fam` containing its subset of samples but all SNPs. -->

---

## 2. Real‐World Datasets

Our “horizontal” assumption holds: each center holds the complete SNP set but has a different subset of participants.

### 2.1. 1000 Genomes Project (1kGP)

- **Samples:** ~2 504  
- **SNPs:** ~chrm5 (phased/imputed VCF).  
- **Centers (Clients):** 5 super‐populations:  
  - AFR, EUR, EAS, SAS, AMR.  
- **Horizontal Splits:**  
  Each super‐population center receives only its population’s samples, but all SNPs (same `.bim` file across centers).  

- **Use Cases:**  
  - Quick functional correctness—end-to-end run.  
  - Demonstrate QC and LR across ancestry‐partitioned centers.

### 2.2. All of Us (AoU)

We simulate realistic “hub/clinic” partitions where each center has different subsets of participants but all SNPs available. 

1. **By Enrollment Site / Health Center**  
   - Partition participants by real clinic/site codes (if existed).  
   - Each center gets all SNPs but only the participants who enrolled at that site.

2. **By Phenotype/Cohort Group**  
   - Partition participants by disease cohort or ancestry group (e.g. diabetes cohort vs. control) as separate centers.  
   - Each center still has the full SNP matrix (after merging local data), but a different set of sample IDs.

3. **Uneven Sample Sizes**  
   - Some centers have large enrollments (e.g. 50 k), others have smaller (e.g. 10 k), reflecting real service‐area populations.  
   - All share the same SNP set.  

- **Clients (Centers):** 7 (simulate distinct data‐collection sites).  
- **Use Cases:**  
  - Demonstrate federated QC and association screening across heterogeneous sample sizes.  
  - Test end‐to‐end runtime on a moderate AoU subset (e.g. 100 k samples, ~10 million SNPs) in ~1–2 days.

---


## 3. Test Scenario Matrix

Below, each scenario specifies (1) dataset scale, (2) partition strategy, (3) number of clients, (4) machines, and (5) expected runtime. All “partition strategy” entries assume horizontal partitioning of SNPs; only sample sets (and in some cases phenotype prevalence) vary.

| Scenario ID | Dataset       | Partition Strategy               | #Clients | Machines            | Expected Runtime   | Notes                                                                             |
|:-----------:|:-------------:|:--------------------------------:|:--------:|:-------------------:|:------------------:|:---------------------------------------------------------------------------------:|
| **A1**      | 1kGP          | Even Sample Split, balanced phenotype | 5        | 1 VM (5 clients)    | ~30 min            | Functional correctness on real‐world data; balanced case/control in each center  |
| **A2**      | Sim‐Tiny      | Even Sample Split, balanced phenotype | 2        | 1 VM (2 clients)    | ~5 min             | Quick sanity correctness                                                           |
| **B1**      | Sim‐Small     | Even Sample Split, balanced phenotype | 3        | 1 VM (3 clients)    | ~1 h               | Small synthetic, equal partitions                                                 |
| **B2**      | Sim‐Small     | Uneven Sample Sizes, balanced phenotype | 3     | 1 VM (3 clients)    | ~1.5 h             | Impact of sample‐size imbalance                                                    |
| **B3**      | Sim‐Small     | Uneven Size + Phenotype Skew      | 3        | 1 VM (3 clients)    | ~1.5 h             | Center A large balanced, B medium moderately skewed, C small highly skewed         |
| **C1**      | Sim‐Medium    | Even Sample Split, balanced phenotype | 5     | 5 VMs (1 client/VM) | ~4 h               | Medium synthetic, distributed                                                    |
| **C2**      | Sim‐Medium    | Uneven Sample Sizes, balanced phenotype | 5   | 5 VMs (1 client/VM) | ~5 h               | Sample‐size skew across centers                                                   |
| **C3**      | Sim‐Medium    | Uneven Size + Phenotype Skew      | 5        | 5 VMs (1 client/VM) | ~5 h               | Center A: most samples balanced; others smaller & increasingly skewed              |
| **D**      | Sim‐Large     | SF-GWAS partition strategy | 7     | 7 VMs (1 client/VM) | ~1–2 days          | Large synthetic—UKB proxy |
| **E1**      | All of Us (sub)| By Enrollment Site (uneven sizes + real phenotype) | 7 | 7 VMs (1 client/VM) | ~12–24 h          | Real AoU subset (~100 k samples, ~10 M SNPs); clinic‐based splits, phenotype varies by clinic |
| **E2**      | All of Us (sub)| By Disease Cohort (skewed phenotype) | 7  | 7 VMs (1 client/VM) | ~12–24 h          | Each center is a disease cohort (100 % cases or 100 % controls per center); tests extreme skew |
| **E3**      | 1kGP          | Simulate Phenotype Skew Across Populations | 5 | 1 VM (5 clients)    | ~30 min            | 1kGP split by super‐population, then induce phenotype skew in one or two centers   |

> **Explanation of E1–E3:**  
> - E1: Centers correspond to real “health‐center” codes in the AoU metadata; sample counts vary, phenotype prevalence comes from real data.  
> - E2: Centers correspond to disease cohorts—one center may be all cases, another all controls, etc., to test pipeline under extreme phenotype skew.  
> - E3: Take 1kGP’s super‐population splits but artificially skew phenotype (e.g., assign 50/50 cases/controls in EUR, 80/20 in AFR, etc.) to test performance under phenotype biases.

---

## 4. Participation & Dropout Patterns

**Chose one of the datasets and partition strategy**:

Because **late joins are not supported**, every client must begin at **Local QC**. If a client opts out or exits at any stage, it will not rejoin:

1. **Opt‐Out During QC**  
   - Set `participation.local_qc = false`.  
   - Client immediately exits; server proceeds with remaining clients.

2. **Opt‐Out During Iterative KING**  
   - After QC & sync & init_chunks, set `participation.iterative_king = false`.  
   - That client finishes preceding stages, then exits. Remaining clients proceed to LR.

3. **Opt‐Out During Iterative LR**  
   - Client sets `participation.iterative_lr = false`.  
   - Client completes all prior stages, then exits before LR. Remaining clients finish LR.

4. **No Late Joins**  
   - Any client that skips QC (or any earlier stage) is excluded from all downstream stages.



## 5. Test Steps (Per Scenario)

For each Scenario (e.g., “B3: 1 VM, 3 clients, Sim‐Small, Uneven Size + Phenotype Skew”):

1. **Provision Environment**  
   - Build/pull Docker image or VM template with PLINK, Python 3.x, Flower, and dependencies.  
   - Launch VMs/containers as specified (“Machines”).

2. **Generate or Download Data**  
   - **Synthetic (Sim)**:  
     ```bash
     python scripts/gen_synthetic.py \
       --n_samples 5000 \
       --n_snps 50000 \
       --centers 3 \
       --uneven-samples \
       --phenotype-skew \
       --skew-ratios 0.5,0.3,0.1 \      # e.g., 50/50, 30/70, 10/90 case/control ratios
       --out_dir /data/sim_small_skew
     ```  
     - `--uneven-samples` indicates non‐equal sample counts.  
     - `--phenotype-skew` and `--skew-ratios` define how cases/controls are distributed per center (e.g., 50/50, 70/30, 90/10).  

   - **1kGP:**  
     - Use DataLoader’s `load_1kgp()` to download and convert to PLINK.  
     - Split the 2 504 samples by super‐population and then apply phenotype skew flags if simulating skew.  

   - **All of Us:**  
     - Obtain 100 k–200 k sample subset and convert via `load_allofus()`.  
     - Partition by enrollment site or cohort with real phenotype distributions (from metadata) or artificially impose skew.

3. **Configure Clients**  
   - Each client’s `config.yaml` points to its local PLINK prefix and sets `participation` flags:  
     ```yaml
     input_data:
       path: "/data/sim_small_skew/center1/prefix"
       type: "bed"
     output:
       intermediate_dir: "/output/center1/"
       log_dir: "/logs/center1/"
     thresholds:
       maf_threshold: 0.01
       missing_threshold: 0.05
       hwe_threshold: 1e-6
       p_threshold: 1e-3
     participation:
       local_qc: true
       global_qc: true
       sync: true
       init_chunks: true
       iterative_king: true   # set to false if client should exit at KING
       local_lr: true
       local_lr_filter_response: true
       init_chunks_lr: true
       iterative_lr: true
     bypass_prereq: {}
     ```

4. **Start Server**  
   ```bash
   cd server/
   python main_server.py


5. Monitor & Collect Metrics

  - Server logs (participants_per_stage, failures_per_stage) record which clients participated or exited at each stage.
  - Client logs (logs/iteration_log.txt) capture timestamps for stage start/end or exit.
  - Resource Usage:
  - Run /usr/bin/time -v python main_client.py for memory and I/O metrics.
  - Run iftop or nload on the server to log network bandwidth.

6. Validate Outputs

  - Global QC: Compare SNP‐exclusion lists to a single‐machine PLINK QC on the merged dataset.
  - KING: Compare federated kinship estimates to a centralized PLINK KING on combined samples.
  - LR: Compare federated association p-values to centralized PLINK logistic regression.

7. Document Results
 - Populate a spreadsheet or JSON log with columns: ScenarioID, ClientID, Stage, StartTime, EndTime, Duration(s), MemoryPeak(MB), BytesSent(Bytes), BytesReceived(Bytes), ExitStage (if any), ExitReason


## Reporting & Analysis

After all scenarios are complete:
  1. Aggregate all logs into a central directory.
  2. Generate summary tables/visualizations:
  - Per‐stage runtime vs. dataset size (line charts).
  - Client dropout analysis (bar charts showing number of clients finishing each stage).
  - Resource usage per stage (memory, CPU, I/O).
  - Network traffic (bytes sent/received per stage).
  3. Compare to SF‐GWAS benchmarks on UKB‐scale tasks.
  4. Prepare a final report (Markdown or Jupyter notebook) summarizing:
  - Correctness validation (QC, KING, LR results match centralized runs).
  - Robustness (clients exiting at specified stages do not break pipeline).
  - Scalability (timing & resource usage from Tiny → Large).
  - Real‐world applicability (1kGP and AoU subset results).

⸻
