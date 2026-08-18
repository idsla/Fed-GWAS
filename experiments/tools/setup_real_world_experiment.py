#!/usr/bin/env python3
"""
Setup Real-World Experiment

Helper script to set up real-world dataset experiments:
- Validate PLINK data files
- Check for phenotype data
- Generate config files from template
- Partition data across clients (optional)
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional, List, Tuple
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def find_plink() -> str:
    """Find PLINK binary."""
    import shutil
    plink = shutil.which("plink") or shutil.which("plink2")
    if not plink:
        raise RuntimeError("PLINK not found in PATH. Please install PLINK.")
    return plink


def validate_plink_files(plink_prefix: str) -> Tuple[bool, List[str]]:
    """
    Validate PLINK files exist and are readable.
    
    Args:
        plink_prefix: Path to PLINK files (without extension)
        
    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []
    base_path = Path(plink_prefix)
    
    # Check required files
    required_files = ['.bed', '.bim', '.fam']
    for ext in required_files:
        file_path = base_path.with_suffix(ext)
        if not file_path.exists():
            errors.append(f"Missing file: {file_path}")
    
    if errors:
        return False, errors
    
    # Try to read with PLINK
    try:
        plink = find_plink()
        result = subprocess.run(
            [plink, '--bfile', str(base_path), '--check'],
            capture_output=True,
            text=True,
            timeout=60
        )
        if result.returncode != 0:
            errors.append(f"PLINK validation failed: {result.stderr}")
            return False, errors
    except subprocess.TimeoutExpired:
        errors.append("PLINK validation timed out")
        return False, errors
    except Exception as e:
        errors.append(f"PLINK validation error: {e}")
        return False, errors
    
    return True, []


def check_phenotypes(plink_prefix: str) -> Tuple[bool, dict]:
    """
    Check if phenotypes are present in .fam file.
    
    Args:
        plink_prefix: Path to PLINK files
        
    Returns:
        Tuple of (has_phenotypes, phenotype_info)
    """
    fam_file = Path(f"{plink_prefix}.fam")
    if not fam_file.exists():
        return False, {'error': 'FAM file not found'}
    
    try:
        fam_df = pd.read_csv(
            fam_file,
            sep='\t',
            header=None,
            names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype']
        )
        
        # Check phenotype column
        phenotypes = fam_df['Phenotype'].values
        unique_phenos = set(phenotypes)
        
        # Valid phenotypes: 1 (control), 2 (case), -9 or 0 (missing)
        valid_phenos = {1, 2, -9, 0}
        invalid_phenos = unique_phenos - valid_phenos
        
        if invalid_phenos:
            return False, {
                'error': f'Invalid phenotype values: {invalid_phenos}',
                'unique_phenotypes': list(unique_phenos)
            }
        
        # Count cases and controls
        case_count = int((phenotypes == 2).sum())
        control_count = int((phenotypes == 1).sum())
        missing_count = int(((phenotypes == -9) | (phenotypes == 0)).sum())
        
        if case_count == 0 and control_count == 0:
            return False, {
                'error': 'No valid phenotypes found (all missing)',
                'missing_count': missing_count
            }
        
        if case_count == 0:
            return False, {
                'error': 'No cases found (all controls or missing)',
                'control_count': control_count,
                'missing_count': missing_count
            }
        
        if control_count == 0:
            return False, {
                'error': 'No controls found (all cases or missing)',
                'case_count': case_count,
                'missing_count': missing_count
            }
        
        case_rate = case_count / (case_count + control_count) if (case_count + control_count) > 0 else 0
        
        return True, {
            'case_count': case_count,
            'control_count': control_count,
            'missing_count': missing_count,
            'case_rate': case_rate,
            'total_samples': len(fam_df)
        }
    
    except Exception as e:
        return False, {'error': f'Failed to read FAM file: {e}'}


def validate_dataset(dataset_name: str, data_dir: Optional[str] = None) -> bool:
    """
    Validate a real-world dataset setup.
    
    Args:
        dataset_name: Name of dataset (e.g., '1000genomes')
        data_dir: Optional data directory path
        
    Returns:
        True if valid, False otherwise
    """
    if data_dir is None:
        data_dir = f"experiments/real_world/{dataset_name}/data"
    
    data_path = Path(data_dir)
    if not data_path.exists():
        logger.error(f"Data directory not found: {data_path}")
        return False
    
    # Check for genotype files
    logger.info(f"Validating dataset: {dataset_name}")
    
    # Look for PLINK files
    plink_files = list(data_path.glob("*.bed"))
    if not plink_files:
        logger.error("No PLINK .bed files found")
        return False
    
    # Validate each PLINK file set
    all_valid = True
    for bed_file in plink_files:
        plink_prefix = str(bed_file.with_suffix(''))
        logger.info(f"Validating: {plink_prefix}")
        
        # Validate PLINK files
        is_valid, errors = validate_plink_files(plink_prefix)
        if not is_valid:
            logger.error(f"PLINK validation failed: {errors}")
            all_valid = False
            continue
        
        # Check phenotypes
        has_phenos, pheno_info = check_phenotypes(plink_prefix)
        if not has_phenos:
            logger.error(f"Phenotype check failed: {pheno_info}")
            all_valid = False
            continue
        
        logger.info(f"✓ Valid: {plink_prefix}")
        logger.info(f"  Cases: {pheno_info['case_count']}, Controls: {pheno_info['control_count']}")
        logger.info(f"  Case rate: {pheno_info['case_rate']:.3f}")
    
    return all_valid


def create_config_template(dataset_name: str, output_dir: Optional[str] = None) -> None:
    """
    Create a config template for a dataset.
    
    Args:
        dataset_name: Name of dataset
        output_dir: Output directory for config files
    """
    if output_dir is None:
        output_dir = f"experiments/real_world/{dataset_name}"
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Create main config.yaml
    config_content = f"""experiment_name: {dataset_name}
experiment_category: real_world
description: 'Real-world dataset: {dataset_name} with simulated phenotypes'

clients:
  config_files:
    0: experiments/real_world/{dataset_name}/configs/center_1/config.yaml
    1: experiments/real_world/{dataset_name}/configs/center_2/config.yaml

data:
  data_dir: experiments/real_world/{dataset_name}/data
  partition_strategy: even

server:
  num_server_rounds: 15
  chunk_size: 100
  min_available_clients: 1
  min_fit_clients: 1

analysis:
  generate_baseline: true
  compare_results: true
  metrics_collection: true
"""
    
    config_file = output_path / "config.yaml"
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    logger.info(f"Created config template: {config_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Setup and validate real-world dataset experiments'
    )
    parser.add_argument(
        '--dataset',
        type=str,
        required=True,
        help='Dataset name (e.g., 1000genomes)'
    )
    parser.add_argument(
        '--data-dir',
        type=str,
        default=None,
        help='Data directory path (default: experiments/real_world/<dataset>/data)'
    )
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Validate dataset setup'
    )
    parser.add_argument(
        '--create-config',
        action='store_true',
        help='Create config template'
    )
    
    args = parser.parse_args()
    
    if args.validate:
        success = validate_dataset(args.dataset, args.data_dir)
        sys.exit(0 if success else 1)
    
    if args.create_config:
        create_config_template(args.dataset)
        logger.info("Config template created. Edit config files as needed.")
    
    if not args.validate and not args.create_config:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
