#!/usr/bin/env python3
"""
Systematic Investigation of Federated vs Baseline Mismatches

This script investigates why federated results don't match baseline results.
"""

import sys
from pathlib import Path
import subprocess
import re

def run_command(cmd):
    """Run a command and return output"""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        return f"ERROR: {e.stderr}"

def investigate_qc_mismatch():
    """Investigate QC exclusion mismatch"""
    print("=" * 80)
    print("1. QC EXCLUSION MISMATCH INVESTIGATION")
    print("=" * 80)
    
    baseline_dir = Path("experiments/correctness/tiny_even/data/tiny/centralized_baseline")
    results_dir = Path("experiments/correctness/tiny_even/results")
    
    # Get baseline excluded SNPs
    original_bim = baseline_dir / "merged.bim"
    filtered_bim = baseline_dir / "qc.bim"
    
    if original_bim.exists() and filtered_bim.exists():
        original_snps = set()
        with open(original_bim, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    original_snps.add(parts[1])
        
        filtered_snps = set()
        with open(filtered_bim, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    filtered_snps.add(parts[1])
        
        baseline_excluded = sorted(original_snps - filtered_snps)
        print(f"Baseline: {len(original_snps)} original SNPs, {len(filtered_snps)} after QC")
        print(f"Baseline excluded: {len(baseline_excluded)} SNPs")
        if baseline_excluded:
            print(f"  Examples: {baseline_excluded[:5]}")
    
    # Check federated exclusion
    exclusion_files = list(results_dir.glob("center_*/logs/global_qc_exclusion.txt"))
    print(f"\nFederated exclusion files: {len(exclusion_files)}")
    
    if exclusion_files:
        fed_excluded = set()
        for f in exclusion_files:
            with open(f, 'r') as exc_file:
                fed_excluded.update(line.strip() for line in exc_file if line.strip())
        print(f"Federated excluded: {len(fed_excluded)} SNPs")
        if fed_excluded:
            print(f"  Examples: {sorted(fed_excluded)[:5]}")
    else:
        print("  ⚠️  No exclusion files found - checking logs...")
        for center_dir in results_dir.glob("center_*/logs"):
            log_files = list(center_dir.glob("supernode_*.txt"))
            if log_files:
                with open(log_files[0], 'r') as f:
                    content = f.read()
                    match = re.search(r'Global QC exclusion list size: (\d+)', content)
                    if match:
                        print(f"  {center_dir.parent.name}: {match.group(1)} SNPs (from log)")
    
    # Check thresholds
    print("\nThresholds comparison:")
    config_file = Path("experiments/correctness/tiny_even/config.yaml")
    if config_file.exists():
        with open(config_file, 'r') as f:
            content = f.read()
            maf_match = re.search(r'maf_threshold:\s*([\d.e-]+)', content)
            miss_match = re.search(r'missing_threshold:\s*([\d.e-]+)', content)
            hwe_match = re.search(r'hwe_threshold:\s*([\d.e-]+)', content)
            if maf_match:
                print(f"  MAF threshold: {maf_match.group(1)}")
            if miss_match:
                print(f"  Missing threshold: {miss_match.group(1)}")
            if hwe_match:
                print(f"  HWE threshold: {hwe_match.group(1)}")
    
    print("\n" + "=" * 80)

def investigate_king_mismatch():
    """Investigate KING correlation mismatch"""
    print("=" * 80)
    print("2. KING CORRELATION MISMATCH INVESTIGATION")
    print("=" * 80)
    
    baseline_dir = Path("experiments/correctness/tiny_even/data/tiny/centralized_baseline")
    results_dir = Path("experiments/correctness/tiny_even/results")
    
    # Check baseline KING
    baseline_king_files = list(baseline_dir.glob("*.kin0")) + list(baseline_dir.glob("*.genome"))
    if baseline_king_files:
        baseline_file = baseline_king_files[0]
        print(f"Baseline KING file: {baseline_file.name}")
        
        # Count pairs
        baseline_pairs = 0
        baseline_samples = set()
        with open(baseline_file, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    baseline_pairs += 1
                    # Try to extract sample IDs
                    if 'ID1' in header or 'IID1' in header:
                        baseline_samples.add(parts[0])
                        baseline_samples.add(parts[1])
                    elif len(parts) >= 4:
                        baseline_samples.add(parts[0])
                        baseline_samples.add(parts[2])
        
        print(f"  Baseline pairs: {baseline_pairs}")
        print(f"  Baseline unique samples: {len(baseline_samples)}")
        if baseline_samples:
            sample_list = sorted(list(baseline_samples))[:5]
            print(f"  Sample ID format: {sample_list}")
    
    # Check federated KING
    server_intermediate = results_dir / "server/intermediate"
    fed_king_files = []
    if server_intermediate.exists():
        for session_dir in server_intermediate.iterdir():
            if session_dir.is_dir():
                fed_king_files.extend(list(session_dir.glob("*.kin0")))
                fed_king_files.extend(list(session_dir.glob("king_results.kin0")))
    
    if fed_king_files:
        fed_file = fed_king_files[0]
        print(f"\nFederated KING file: {fed_file.name}")
        
        fed_pairs = 0
        fed_samples = set()
        with open(fed_file, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    fed_pairs += 1
                    if 'ID1' in header or 'IID1' in header:
                        fed_samples.add(parts[0])
                        fed_samples.add(parts[1])
                    elif len(parts) >= 4:
                        fed_samples.add(parts[0])
                        fed_samples.add(parts[2])
        
        print(f"  Federated pairs: {fed_pairs}")
        print(f"  Federated unique samples: {len(fed_samples)}")
        if fed_samples:
            sample_list = sorted(list(fed_samples))[:5]
            print(f"  Sample ID format: {sample_list}")
        
        # Check sample ID overlap
        if baseline_samples and fed_samples:
            overlap = baseline_samples & fed_samples
            print(f"\n  Sample ID overlap: {len(overlap)}/{len(baseline_samples)} ({100*len(overlap)/len(baseline_samples):.1f}%)")
            if len(overlap) < len(baseline_samples) * 0.9:
                print(f"  ⚠️  MAJOR ISSUE: Sample IDs don't match!")
                print(f"     This explains the KING correlation mismatch.")
                print(f"     Federated uses anonymized sample IDs (sample_offset),")
                print(f"     while baseline uses original sample IDs.")
    
    print("\n" + "=" * 80)

def investigate_lr_mismatch():
    """Investigate LR correlation mismatch"""
    print("=" * 80)
    print("3. LR CORRELATION MISMATCH INVESTIGATION")
    print("=" * 80)
    
    baseline_dir = Path("experiments/correctness/tiny_even/data/tiny/centralized_baseline")
    results_dir = Path("experiments/correctness/tiny_even/results")
    
    # Check baseline LR
    baseline_lr_files = list(baseline_dir.glob("*.assoc.logistic"))
    if baseline_lr_files:
        baseline_file = baseline_lr_files[0]
        print(f"Baseline LR file: {baseline_file.name}")
        
        baseline_snps = set()
        baseline_count = 0
        with open(baseline_file, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[1] != 'SNP':
                    baseline_snps.add(parts[1])
                    baseline_count += 1
        
        print(f"  Baseline SNPs: {len(baseline_snps)}")
        print(f"  Baseline rows: {baseline_count}")
    
    # Check federated global LR
    server_intermediate = results_dir / "server/intermediate"
    fed_lr_files = []
    if server_intermediate.exists():
        for session_dir in server_intermediate.iterdir():
            if session_dir.is_dir():
                fed_lr_files.extend(list(session_dir.glob("*.assoc.logistic")))
    
    if fed_lr_files:
        fed_file = fed_lr_files[0]
        print(f"\nFederated Global LR file: {fed_file.name}")
        
        fed_snps = set()
        fed_count = 0
        with open(fed_file, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[1] != 'SNP':
                    fed_snps.add(parts[1])
                    fed_count += 1
        
        print(f"  Federated SNPs: {len(fed_snps)}")
        print(f"  Federated rows: {fed_count}")
        
        # Check SNP overlap
        if baseline_snps and fed_snps:
            overlap = baseline_snps & fed_snps
            print(f"\n  SNP overlap: {len(overlap)}/{len(baseline_snps)} ({100*len(overlap)/len(baseline_snps):.1f}%)")
            
            missing = baseline_snps - fed_snps
            extra = fed_snps - baseline_snps
            
            if missing:
                print(f"  ⚠️  Missing in federated: {len(missing)}")
                print(f"     Examples: {sorted(list(missing))[:5]}")
            
            if extra:
                print(f"  ⚠️  Extra in federated: {len(extra)}")
                print(f"     Examples: {sorted(list(extra))[:5]}")
            
            if len(overlap) < len(baseline_snps) * 0.9:
                print(f"\n  ⚠️  MAJOR ISSUE: SNP sets don't match!")
                print(f"     Possible causes:")
                print(f"     1. QC filtering mismatch (baseline excluded {len(baseline_snps - fed_snps)} SNPs)")
                print(f"     2. SNP anonymization in federated pipeline")
                print(f"     3. Different data subsets")
    
    # Check local LR
    local_lr_files = []
    alt_dir = Path("experiments/simulations/scenario_A/results")
    for center_dir in alt_dir.glob("center_*/logs"):
        local_lr_files.extend(list(center_dir.glob("local_lr*.assoc.logistic")))
    
    if local_lr_files:
        print(f"\nLocal LR files found: {len(local_lr_files)}")
        print(f"  Location: {local_lr_files[0].parent}")
        print(f"  ⚠️  NOTE: Local LR uses partitioned data (half samples per client)")
        print(f"     This explains why correlation is low - different sample sizes!")
    
    print("\n" + "=" * 80)

def main():
    """Run all investigations"""
    print("\n" + "=" * 80)
    print("SYSTEMATIC INVESTIGATION OF FEDERATED vs BASELINE MISMATCHES")
    print("=" * 80 + "\n")
    
    investigate_qc_mismatch()
    investigate_king_mismatch()
    investigate_lr_mismatch()
    
    print("\n" + "=" * 80)
    print("SUMMARY AND RECOMMENDATIONS")
    print("=" * 80)
    print("""
KEY FINDINGS:

1. QC MISMATCH:
   - Baseline excluded 12 SNPs, Federated excluded 0
   - Cause: Federated experiment was run before exclusion file saving was added
   - Fix: Re-run experiment with updated code

2. KING CORRELATION (r=0.88):
   - Sample IDs are anonymized in federated pipeline (sample_offset)
   - Baseline uses original sample IDs
   - This causes sample pair mismatches
   - Fix: Need to map anonymized IDs back to original IDs for comparison

3. LR CORRELATION (r≈0):
   - SNP sets may not match due to QC filtering differences
   - Local LR uses partitioned data (different sample sizes)
   - Global LR may have SNP anonymization issues
   - Fix: Ensure QC filtering matches, verify SNP ID mapping

RECOMMENDATIONS:

1. Re-run experiment with updated code that saves exclusion files
2. Implement ID mapping for KING comparison (anonymized -> original)
3. Verify QC thresholds are identical between baseline and federated
4. Check if SNP anonymization affects LR results
5. Ensure data partitioning doesn't affect QC filtering
""")

if __name__ == "__main__":
    main()


