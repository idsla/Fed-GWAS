#!/usr/bin/env python3
"""
Centralized Baseline Generator for Federated GWAS Experiments

This module generates centralized PLINK baselines for comparison:
- Merge federated datasets from all centers
- Run centralized PLINK QC, KING, and LR
- Store results for comparison with federated results
"""

import os
import subprocess
import logging
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any
import numpy as np
import yaml

class BaselineGenerator:
    """Centralized baseline generator class"""
    
    def __init__(self, experiment_config: str, output_dir: Optional[str] = None):
        """Initialize baseline generator"""
        self.config_file = Path(experiment_config)
        self.config = self._load_config()
        
        # Determine output directory
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            # Use experiment data directory
            data_dir = Path(self.config.get('data', {}).get('data_dir', 'experiments/simulations/scenario_A/data/tiny'))
            self.output_dir = data_dir / "centralized_baseline"
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger = self._setup_logger()
        
        # PLINK binary path
        self.plink_binary = self._find_plink()
        
        self.logger.info(f"Baseline generator initialized")
        self.logger.info(f"  Config: {experiment_config}")
        self.logger.info(f"  Output: {self.output_dir}")
    
    def _load_config(self) -> Dict[str, Any]:
        """Load experiment configuration"""
        with open(self.config_file, 'r') as f:
            config = yaml.safe_load(f) or {}
        data = config.get("data")
        if isinstance(data, dict) and data.get("data_dir"):
            data_dir = Path(str(data["data_dir"]))
            if not data_dir.is_absolute():
                data["data_dir"] = str((self.config_file.parent / data_dir).resolve())
        return config

    def _thresholds(self) -> Dict[str, Any]:
        """Return baseline thresholds with top-level values overriding client defaults."""
        thresholds: Dict[str, Any] = {}
        clients = self.config.get("clients", {})
        if isinstance(clients, dict) and isinstance(clients.get("thresholds"), dict):
            thresholds.update(clients["thresholds"])
        if isinstance(self.config.get("thresholds"), dict):
            thresholds.update(self.config["thresholds"])
        return thresholds
    
    def _setup_logger(self) -> logging.Logger:
        """Set up logging"""
        logger = logging.getLogger('baseline_generator')
        logger.setLevel(logging.INFO)
        
        # File handler
        log_file = self.output_dir / "baseline_generation.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        logger.addHandler(file_handler)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        logger.addHandler(console_handler)
        
        return logger
    
    def _find_plink(self) -> str:
        """Find PLINK binary"""
        import platform
        
        project_root = Path(__file__).parent.parent.parent
        
        # Determine OS-specific PLINK path
        system = platform.system().lower()
        if system == "darwin":  # macOS
            plink_dir = project_root / "plink" / "plink_mac"
            plink_name = "plink"
        elif system == "linux":
            plink_dir = project_root / "plink" / "plink_linux"
            plink_name = "plink"
        elif system == "windows":
            plink_dir = project_root / "plink" / "plink_win"
            plink_name = "plink.exe"
        else:
            plink_dir = project_root / "plink" / "plink_mac"  # Default to macOS
            plink_name = "plink"
        
        # Check OS-specific location first
        plink_path = plink_dir / plink_name
        if plink_path.exists() and os.access(plink_path, os.X_OK):
            return str(plink_path)
        
        # Check project bin directory
        bin_plink = project_root / "bin" / plink_name
        if bin_plink.exists() and os.access(bin_plink, os.X_OK):
            return str(bin_plink)
        
        # Try in PATH
        if shutil.which("plink"):
            return "plink"
        if shutil.which("plink2"):
            return "plink2"
        
        # Last resort: return the expected path (will fail with clear error)
        self.logger.warning(f"PLINK not found in expected location: {plink_path}")
        return str(plink_path)

    def _run_plink(self, cmd: List[str], error_msg: str) -> None:
        """Run a PLINK command and raise with stderr on failure."""
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            stderr = e.stderr or e.stdout or str(e)
            self.logger.error(f"{error_msg}: {stderr}")
            raise

    def _apply_mind(self, input_prefix: str, output_prefix: str, mind_threshold: float) -> str:
        """Apply sample-level missingness filter (per-center, matches federated local_qc)."""
        cmd = [
            self.plink_binary,
            "--bfile", input_prefix,
            "--mind", str(mind_threshold),
            "--make-bed",
            "--allow-no-sex",
            "--out", output_prefix
        ]
        self._run_plink(cmd, f"Sample QC (--mind) failed for {input_prefix}")
        return output_prefix
    
    def merge_datasets(self) -> str:
        """Merge all center datasets into a single centralized dataset"""
        self.logger.info("Merging federated datasets...")
        
        data_dir = Path(self.config.get('data', {}).get('data_dir', ''))
        if not data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {data_dir}")
        
        # Find all center directories (ignore sample list files)
        center_dirs = sorted([p for p in data_dir.glob("center_*") if p.is_dir()])
        if not center_dirs:
            raise ValueError(f"No center directories found in {data_dir}")
        
        self.logger.info(f"Found {len(center_dirs)} center directories")
        
        # Apply per-center --mind before merge (matches federated local_qc)
        thresholds = self._thresholds()
        mind_threshold = thresholds.get('local_mind_threshold', thresholds.get('mind_threshold', 0.1))
        mind_root = self.output_dir / "mind_centers"
        mind_root.mkdir(parents=True, exist_ok=True)

        # Get PLINK files from first center
        first_center = center_dirs[0]
        center_files = list(first_center.glob("*.bed"))
        if not center_files:
            raise ValueError(f"No .bed files found in {first_center}")
        
        first_prefix = center_files[0].stem.replace("_center_1", "").replace("_center_2", "").replace("_center_3", "")
        first_bed = first_center / f"{first_prefix}_center_{first_center.name.split('_')[1]}.bed"
        
        if not first_bed.exists():
            # Try alternative naming
            bed_files = list(first_center.glob("*.bed"))
            if bed_files:
                first_bed = bed_files[0]
            else:
                raise ValueError(f"Could not find bed file in {first_center}")
        
        first_prefix = first_bed.stem

        # Apply mind to first center
        first_mind_prefix = mind_root / f"{first_center.name}_mind"
        self.logger.info(f"Applying --mind to {first_center.name} (threshold={mind_threshold})")
        self._apply_mind(str(first_bed.with_suffix('')), str(first_mind_prefix), mind_threshold)
        
        # Copy first dataset as starting point
        merged_prefix = self.output_dir / "merged"
        shutil.copy(f"{first_mind_prefix}.bed", f"{merged_prefix}.bed")
        shutil.copy(f"{first_mind_prefix}.bim", f"{merged_prefix}.bim")
        shutil.copy(f"{first_mind_prefix}.fam", f"{merged_prefix}.fam")
        
        self.logger.info(f"Starting merge with {first_prefix}")
        
        # Merge remaining centers
        for center_dir in center_dirs[1:]:
            center_num = center_dir.name.split('_')[1]
            # Find any bed file in this center
            bed_files = list(center_dir.glob("*.bed"))
            if not bed_files:
                self.logger.warning(f"Skipping {center_dir}: no bed file found")
                continue
            center_bed = bed_files[0]
            center_prefix = center_bed.stem
            
            # Merge using PLINK with multi-stage error recovery
            center_prefix_str = str(center_bed.with_suffix(''))
            center_mind_prefix = mind_root / f"{center_dir.name}_mind"
            self.logger.info(f"Applying --mind to {center_dir.name} (threshold={mind_threshold})")
            self._apply_mind(center_prefix_str, str(center_mind_prefix), mind_threshold)
            merge_tmp = str(self.output_dir / "merge_tmp")
            
            self.logger.info(f"Merging {center_prefix}...")
            
            # Strict merge only (no strand flips or ID disambiguation)
            merge_cmd = [
                self.plink_binary,
                "--bfile", str(merged_prefix),
                "--bmerge", str(center_mind_prefix),
                "--make-bed",
                "--allow-no-sex",
                "--out", merge_tmp
            ]
            
            self._run_plink(merge_cmd, f"Strict merge failed for {center_prefix} (no strand/ID fixes allowed)")
            
            # Move merged files if merge succeeded
            shutil.move(f"{merge_tmp}.bed", f"{merged_prefix}.bed")
            shutil.move(f"{merge_tmp}.bim", f"{merged_prefix}.bim")
            shutil.move(f"{merge_tmp}.fam", f"{merged_prefix}.fam")

            # Clean up temporary merge log
            tmp_log = self.output_dir / "merge_tmp.log"
            if tmp_log.exists():
                tmp_log.unlink()
            
            self.logger.info(f"Successfully merged {center_prefix}")
        
        self.logger.info(f"Merged dataset created: {merged_prefix}")
        return str(merged_prefix)
    
    def run_qc(self, input_prefix: str) -> None:
        """Run centralized QC analysis"""
        self.logger.info("Running centralized QC...")
        
        thresholds = self._thresholds()
        maf_threshold = thresholds.get('maf_threshold', 0.01)
        missing_threshold = thresholds.get('missing_threshold', 0.05)
        hwe_threshold = thresholds.get('hwe_threshold', 1e-6)
        qc_prefix = self.output_dir / "qc"
        qc_input_prefix = str(input_prefix)

        # Compute exclusion list using the same chi-square logic as federated QC
        from pipeline.clients.local_qc import compute_genotype_counts, compute_missingness_counts
        from pipeline.clients.client_qc_aggregator import _compute_exclusion_list

        counts = compute_genotype_counts(
            qc_input_prefix,
            "baseline",
            log_dir=str(self.output_dir),
            plink_binary=self.plink_binary,
        )
        missing = compute_missingness_counts(
            qc_input_prefix,
            "baseline",
            log_dir=str(self.output_dir),
            plink_binary=self.plink_binary,
        )
        thresholds_arr = np.array([maf_threshold, missing_threshold, hwe_threshold], dtype=np.float64)
        exclusion_indices = sorted(list(_compute_exclusion_list(counts, missing, thresholds_arr)))

        # Map indices -> SNP IDs
        bim_file = f"{qc_input_prefix}.bim"
        snp_ids: List[str] = []
        with open(bim_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    snp_ids.append(parts[1])
        excluded_snp_ids = [snp_ids[i] for i in exclusion_indices if 0 <= i < len(snp_ids)]

        # Write exclusion list and apply it
        exclude_file = self.output_dir / "qc_exclude_snps.txt"
        with open(exclude_file, "w") as f:
            for snp_id in excluded_snp_ids:
                f.write(f"{snp_id}\n")

        qc_cmd = [
            self.plink_binary,
            "--bfile", qc_input_prefix,
            "--exclude", str(exclude_file),
            "--make-bed",
            "--allow-no-sex",
            "--out", str(qc_prefix)
        ]
        self._run_plink(qc_cmd, "QC exclusion (chi-square) failed")
        self.logger.info("QC analysis completed (chi-square exclusion list applied)")
        
        # Also run missingness check
        missing_cmd = [
            self.plink_binary,
            "--bfile", str(qc_prefix),
            "--missing",
            "--out", str(qc_prefix)
        ]
        try:
            subprocess.run(missing_cmd, check=True, capture_output=True, text=True)
            self.logger.info("Missingness analysis completed")
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"Missingness check failed: {e}")
        
        # Run frequency analysis for MAF
        freq_cmd = [
            self.plink_binary,
            "--bfile", str(qc_prefix),
            "--freq",
            "--out", str(qc_prefix)
        ]
        try:
            subprocess.run(freq_cmd, check=True, capture_output=True, text=True)
            self.logger.info("Frequency analysis completed")
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"Frequency check failed: {e}")
        
        # Run HWE test
        hwe_cmd = [
            self.plink_binary,
            "--bfile", str(qc_prefix),
            "--hardy",
            "--out", str(qc_prefix)
        ]
        try:
            subprocess.run(hwe_cmd, check=True, capture_output=True, text=True)
            self.logger.info("HWE test completed")
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"HWE test failed: {e}")
    
    def run_king(self, input_prefix: str) -> None:
        """Run centralized KING kinship estimation"""
        self.logger.info("Running centralized KING...")
        
        king_prefix = self.output_dir / "king"
        
        # Try PLINK2 first (preferred for KING)
        import platform
        system = platform.system().lower()
        project_root = Path(__file__).parent.parent.parent
        
        # Find plink2 in same directory as plink
        if system == "darwin":
            plink2_dir = project_root / "plink" / "plink_mac"
        elif system == "linux":
            plink2_dir = project_root / "plink" / "plink_linux"
        elif system == "windows":
            plink2_dir = project_root / "plink" / "plink_win"
            plink2_name = "plink2.exe"
        else:
            plink2_dir = project_root / "plink" / "plink_mac"
        
        if system != "windows":
            plink2_name = "plink2"
        
        plink2_binary = plink2_dir / plink2_name
        
        # Check if plink2 exists
        if plink2_binary.exists() and os.access(plink2_binary, os.X_OK):
            plink2_binary = str(plink2_binary)
        elif shutil.which("plink2"):
            plink2_binary = "plink2"
        else:
            plink2_binary = None
        
        if plink2_binary:
            king_cmd = [
                plink2_binary,
                "--bfile", input_prefix,
                "--make-king-table",
                "--out", str(king_prefix)
            ]
        else:
            # Fall back to PLINK 1.9 genome command
            self.logger.warning("PLINK2 not found, using PLINK 1.9 --genome")
            king_cmd = [
                self.plink_binary,
                "--bfile", input_prefix,
                "--genome",
                "--out", str(king_prefix)
            ]
        
        try:
            subprocess.run(king_cmd, check=True, capture_output=True, text=True)
            self.logger.info("KING analysis completed")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"KING failed: {e}")
            raise

    def filter_samples_by_king(self, input_prefix: str) -> Optional[str]:
        """Filter related samples based on centralized KING results.

        Uses the same threshold as federated clients (thresholds.king_threshold).
        Returns the filtered prefix if samples were removed, otherwise None.
        """
        thresholds = self.config.get('thresholds', {})
        king_threshold = thresholds.get('king_threshold', 0.4)
        try:
            king_threshold = float(king_threshold)
        except Exception:
            king_threshold = 0.4

        king_file = self.output_dir / "king.kin0"
        if not king_file.exists():
            self.logger.warning("KING result file not found; skipping KING-based filtering")
            return None

        # Load related pairs from KING results
        related_samples = set()
        with king_file.open("r") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) < 8:
                    continue
                fid1, fid2 = parts[0], parts[2]
                try:
                    kin = float(parts[7])
                except ValueError:
                    continue
                if kin > king_threshold:
                    related_samples.add(fid1)
                    related_samples.add(fid2)

        if not related_samples:
            self.logger.info(f"No samples exceed KING threshold {king_threshold}; skipping KING-based filtering")
            return None

        # Build FID/IID mapping from input .fam
        fam_file = f"{input_prefix}.fam"
        fid_to_iid = {}
        try:
            with open(fam_file, "r") as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        fid, iid = parts[0], parts[1]
                        fid_to_iid[fid] = iid
        except FileNotFoundError:
            self.logger.warning(f"Cannot find FAM file for KING filtering: {fam_file}")
            return None

        remove_file = self.output_dir / "king_remove_samples.txt"
        removed = 0
        with remove_file.open("w") as f:
            for fid in sorted(related_samples):
                iid = fid_to_iid.get(fid, fid)
                f.write(f"{fid}\t{iid}\n")
                removed += 1

        if removed == 0:
            self.logger.info("KING remove list empty; skipping KING-based filtering")
            return None

        filtered_prefix = str(self.output_dir / "king_filtered")
        cmd = [
            self.plink_binary,
            "--bfile", input_prefix,
            "--remove", str(remove_file),
            "--make-bed",
            "--allow-no-sex",
            "--out", filtered_prefix,
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info(
                f"KING filtering complete: removed {removed} samples with kinship > {king_threshold}"
            )
            return filtered_prefix
        except subprocess.CalledProcessError as e:
            self.logger.error(f"KING-based filtering failed: {e}")
            return None
    
    def run_lr(self, input_prefix: str) -> None:
        """Run centralized logistic regression"""
        self.logger.info("Running centralized logistic regression...")
        
        lr_prefix = self.output_dir / "lr"
        
        lr_cmd = [
            self.plink_binary,
            "--bfile", input_prefix,
            "--logistic",
            "--allow-no-sex",  # Allow analysis when sex information is missing/ambiguous
            "--out", str(lr_prefix)
        ]
        
        try:
            subprocess.run(lr_cmd, check=True, capture_output=True, text=True)
            self.logger.info("Logistic regression completed")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"LR failed: {e}")
            # Check if it's a phenotype issue
            if "phenotype" in str(e).lower() or "constant" in str(e).lower():
                self.logger.warning("Phenotype may be constant. Try using --pheno or checking phenotype file.")
            raise
    
    def generate_baseline(self) -> Dict[str, str]:
        """Generate complete centralized baseline"""
        self.logger.info("Starting baseline generation...")
        
        results = {}
        
        # Step 1: Merge datasets
        merged_prefix = self.merge_datasets()
        results['merged_prefix'] = merged_prefix
        
        # Step 2: Run QC
        self.run_qc(merged_prefix)
        qc_prefix = self.output_dir / "qc"
        results['qc_prefix'] = str(qc_prefix)
        
        # Step 3: Run KING (prefer QC-filtered data if available)
        if qc_prefix.with_suffix('.bed').exists():
            self.logger.info("Running KING on QC-filtered data")
            self.run_king(str(qc_prefix))
        else:
            self.logger.info("QC data not found; running KING on merged data")
            self.run_king(merged_prefix)
        results['king_prefix'] = str(self.output_dir / "king")
        
        # Step 3b: Optional KING-based sample filtering (for LR consistency)
        king_filtered_prefix = None
        if qc_prefix.with_suffix('.bed').exists():
            king_filtered_prefix = self.filter_samples_by_king(str(qc_prefix))
        else:
            king_filtered_prefix = self.filter_samples_by_king(merged_prefix)
        if king_filtered_prefix:
            results['king_filtered_prefix'] = king_filtered_prefix
        
        # Step 4: Run LR (on QC-filtered data)
        if king_filtered_prefix:
            self.run_lr(king_filtered_prefix)
            results['lr_prefix'] = str(self.output_dir / "lr")
        elif qc_prefix.with_suffix('.bed').exists():
            self.run_lr(str(qc_prefix))
            results['lr_prefix'] = str(self.output_dir / "lr")
        else:
            # Fall back to merged data if QC didn't create new files
            self.run_lr(merged_prefix)
            results['lr_prefix'] = str(self.output_dir / "lr")
        
        self.logger.info("Baseline generation completed!")
        self.logger.info(f"Results saved to: {self.output_dir}")
        
        return results


def main():
    """Command-line interface for baseline generation"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate centralized PLINK baseline for comparison")
    parser.add_argument("experiment_config", type=str, 
                       help="Path to experiment config YAML file")
    parser.add_argument("--output", type=str, 
                       help="Output directory for baseline results (default: data_dir/centralized_baseline)")
    
    args = parser.parse_args()
    
    generator = BaselineGenerator(args.experiment_config, args.output)
    results = generator.generate_baseline()
    
    print("\nBaseline generation completed!")
    print(f"\nResults:")
    for key, value in results.items():
        print(f"  {key}: {value}")


def generate_baseline(experiment_config: str | Path, output_dir: str | Path | None = None) -> Dict[str, str]:
    """Generate a centralized baseline through the core pipeline tool.

    Args:
        experiment_config: Path to a FedGWAS experiment `config.yaml`.
        output_dir: Optional output directory for baseline artifacts.

    Returns:
        Mapping of artifact labels to generated PLINK prefixes.
    """
    generator = BaselineGenerator(
        str(experiment_config),
        str(output_dir) if output_dir is not None else None,
    )
    return generator.generate_baseline()


if __name__ == "__main__":
    main()
