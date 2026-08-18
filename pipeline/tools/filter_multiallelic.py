#!/usr/bin/env python3
"""
Filter Multiallelic Variants from PLINK Data

Standard GWAS practice: Remove variants with 3+ alleles before analysis.
PLINK 1.9 does not handle multiallelic variants well, so they should be excluded
early in the preprocessing pipeline.

This script uses PLINK2's --max-alleles flag (if available) or PLINK 1.9's
--freq command to detect and exclude multiallelic variants.

Usage:
    python filter_multiallelic.py \
        --plink-prefix data/genotypes \
        --output-prefix data/genotypes_biallelic \
        [--keep-log]
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path
import shutil

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def find_plink_binary():
    """Find PLINK binary."""
    import platform
    
    # Try PLINK2 first (has --max-alleles flag)
    plink2 = shutil.which("plink2")
    if plink2:
        return plink2, "plink2"
    
    # Try PLINK 1.9
    plink = shutil.which("plink")
    if plink:
        return plink, "plink"
    
    # Try to find in project plink directory
    project_root = Path(__file__).parent.parent.parent
    system = platform.system().lower()
    
    # Try PLINK2 first
    if system == "darwin":
        plink2_path = project_root / "plink" / "plink2_mac" / "plink2"
        plink_path = project_root / "plink" / "plink_mac" / "plink"
    elif system == "linux":
        plink2_path = project_root / "plink" / "plink2_linux" / "plink2"
        plink_path = project_root / "plink" / "plink_linux" / "plink"
    elif system == "windows":
        plink2_path = project_root / "plink" / "plink2_win" / "plink2.exe"
        plink_path = project_root / "plink" / "plink_win" / "plink.exe"
    else:
        plink2_path = project_root / "plink" / "plink2_mac" / "plink2"
        plink_path = project_root / "plink" / "plink_mac" / "plink"
    
    if plink2_path.exists():
        return str(plink2_path), "plink2"
    if plink_path.exists():
        return str(plink_path), "plink"
    
    raise RuntimeError("PLINK not found. Please install PLINK or ensure it's in PATH.")


def filter_multiallelic_plink2(plink_prefix: str, output_prefix: str) -> bool:
    """
    Use PLINK2's --max-alleles flag to filter multiallelic variants.
    
    This is the most direct and reliable method.
    
    Args:
        plink_prefix: Input PLINK prefix
        output_prefix: Output PLINK prefix (filtered)
        
    Returns:
        True if successful, False otherwise
    """
    plink_binary, plink_type = find_plink_binary()
    
    if plink_type != "plink2":
        return False
    
    logger.info("Using PLINK2 --max-alleles 2 to filter multiallelic variants")
    
    cmd = [
        plink_binary,
        "--bfile", plink_prefix,
        "--max-alleles", "2",
        "--make-bed",
        "--allow-no-sex",
        "--out", output_prefix
    ]
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        # Parse output to count excluded variants
        if result.stdout:
            for line in result.stdout.split('\n'):
                if 'variants removed' in line.lower() or 'variants excluded' in line.lower():
                    logger.info(f"  {line.strip()}")
        
        logger.info(f"✓ Successfully filtered multiallelic variants using PLINK2")
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"PLINK2 command failed: {e}")
        if e.stderr:
            logger.error(f"STDERR: {e.stderr[:500]}")
        return False


def filter_multiallelic_plink19(plink_prefix: str, output_prefix: str) -> bool:
    """
    Use PLINK 1.9 to filter multiallelic variants.
    
    PLINK 1.9 doesn't have a direct --max-alleles flag, so we use a workaround:
    1. Run --freq which will report errors for multiallelic variants
    2. Parse the log to identify problematic variants
    3. Exclude them using --exclude
    
    However, a simpler approach is to use PLINK's --vcf conversion which
    automatically handles multiallelic variants. But since we already have
    PLINK format, we'll use --freq + error parsing.
    
    Args:
        plink_prefix: Input PLINK prefix
        output_prefix: Output PLINK prefix (filtered)
        
    Returns:
        True if successful, False otherwise
    """
    plink_binary, plink_type = find_plink_binary()
    
    logger.info("Using PLINK 1.9 --freq to detect multiallelic variants")
    
    # Step 1: Run --freq to detect multiallelic variants
    freq_prefix = f"{output_prefix}_freq_check"
    cmd_freq = [
        plink_binary,
        "--bfile", plink_prefix,
        "--freq",
        "--allow-no-sex",
        "--out", freq_prefix
    ]
    
    multiallelic_snps = []
    
    try:
        result = subprocess.run(
            cmd_freq,
            capture_output=True,
            text=True
        )
        
        # Parse log file for multiallelic errors
        log_file = Path(f"{freq_prefix}.log")
        if log_file.exists():
            with open(log_file, 'r') as f:
                for line in f:
                    # PLINK reports multiallelic variants in error messages
                    if '3+ alleles' in line.lower() or 'multiallelic' in line.lower():
                        # Try to extract SNP ID from the line
                        parts = line.split()
                        for part in parts:
                            # SNP IDs are typically alphanumeric
                            if len(part) > 2 and not part.startswith('chr') and ':' not in part:
                                multiallelic_snps.append(part)
        
        # Also check stderr for multiallelic errors
        if result.stderr:
            for line in result.stderr.split('\n'):
                if '3+ alleles' in line.lower() or 'multiallelic' in line.lower():
                    parts = line.split()
                    for part in parts:
                        if len(part) > 2 and not part.startswith('chr') and ':' not in part:
                            multiallelic_snps.append(part)
        
        # Clean up freq files
        for ext in ['.freq', '.log', '.nosex']:
            freq_file = Path(f"{freq_prefix}{ext}")
            if freq_file.exists():
                freq_file.unlink()
        
    except Exception as e:
        logger.warning(f"Error during --freq check: {e}")
        # Continue anyway - may not have multiallelic variants
    
    # Step 2: If multiallelic variants detected, exclude them
    if multiallelic_snps:
        # Remove duplicates
        multiallelic_snps = list(set(multiallelic_snps))
        logger.info(f"Detected {len(multiallelic_snps)} multiallelic variants")
        
        exclude_file = Path(f"{output_prefix}_multiallelic_exclude.txt")
        with open(exclude_file, 'w') as f:
            for snp in multiallelic_snps:
                f.write(f"{snp}\n")
        
        cmd_exclude = [
            plink_binary,
            "--bfile", plink_prefix,
            "--exclude", str(exclude_file),
            "--make-bed",
            "--allow-no-sex",
            "--quiet",
            "--out", output_prefix
        ]
        
        try:
            subprocess.run(cmd_exclude, check=True, capture_output=True, text=True)
            exclude_file.unlink()
            logger.info(f"✓ Excluded {len(multiallelic_snps)} multiallelic variants")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to exclude multiallelic variants: {e}")
            if exclude_file.exists():
                exclude_file.unlink()
            return False
    else:
        # No multiallelic variants detected, just copy files
        logger.info("No multiallelic variants detected, copying files as-is")
        for ext in ['.bed', '.bim', '.fam']:
            src = Path(f"{plink_prefix}{ext}")
            dst = Path(f"{output_prefix}{ext}")
            if src.exists():
                shutil.copy2(src, dst)
        return True


def filter_multiallelic_variants(plink_prefix: str, output_prefix: str, keep_log: bool = False) -> bool:
    """
    Main function to filter multiallelic variants from PLINK data.
    
    Standard GWAS practice:
    1. Use PLINK2 --max-alleles 2 if available (most reliable)
    2. Otherwise use PLINK 1.9 --freq to detect and exclude
    
    Args:
        plink_prefix: Input PLINK prefix (path without .bed/.bim/.fam)
        output_prefix: Output PLINK prefix for filtered data
        keep_log: Whether to keep PLINK log file
        
    Returns:
        True if successful, False otherwise
    """
    logger.info("=" * 80)
    logger.info("Filtering Multiallelic Variants (Standard GWAS Preprocessing)")
    logger.info("=" * 80)
    logger.info(f"Input:  {plink_prefix}")
    logger.info(f"Output: {output_prefix}")
    
    # Try PLINK2 first (most reliable method)
    plink_binary, plink_type = find_plink_binary()
    
    if plink_type == "plink2":
        success = filter_multiallelic_plink2(plink_prefix, output_prefix)
        if success:
            # Count variants before and after
            input_bim = Path(f"{plink_prefix}.bim")
            output_bim = Path(f"{output_prefix}.bim")
            
            if input_bim.exists() and output_bim.exists():
                input_count = sum(1 for _ in open(input_bim))
                output_count = sum(1 for _ in open(output_bim))
                excluded_count = input_count - output_count
                
                logger.info("=" * 80)
                logger.info("Filtering Summary:")
                logger.info(f"  Input variants:  {input_count:,}")
                logger.info(f"  Excluded:        {excluded_count:,}")
                logger.info(f"  Output variants: {output_count:,}")
                logger.info("=" * 80)
            
            # Clean up log unless requested to keep
            if not keep_log:
                log_file = Path(f"{output_prefix}.log")
                if log_file.exists():
                    log_file.unlink()
            
            return True
        else:
            logger.warning("PLINK2 filtering failed, falling back to PLINK 1.9 method")
    
    # Fall back to PLINK 1.9 method
    success = filter_multiallelic_plink19(plink_prefix, output_prefix)
    
    if success:
        # Count variants before and after
        input_bim = Path(f"{plink_prefix}.bim")
        output_bim = Path(f"{output_prefix}.bim")
        
        if input_bim.exists() and output_bim.exists():
            input_count = sum(1 for _ in open(input_bim))
            output_count = sum(1 for _ in open(output_bim))
            excluded_count = input_count - output_count
            
            logger.info("=" * 80)
            logger.info("Filtering Summary:")
            logger.info(f"  Input variants:  {input_count:,}")
            logger.info(f"  Excluded:        {excluded_count:,}")
            logger.info(f"  Output variants: {output_count:,}")
            logger.info("=" * 80)
        
        # Clean up log unless requested to keep
        if not keep_log:
            log_file = Path(f"{output_prefix}.log")
            if log_file.exists():
                log_file.unlink()
    
    return success


def main():
    parser = argparse.ArgumentParser(
        description='Filter multiallelic variants from PLINK data (standard GWAS preprocessing)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Standard GWAS Practice:
  Multiallelic variants (3+ alleles) should be excluded before analysis because:
  1. PLINK 1.9 does not handle them well during merging
  2. Most GWAS methods assume biallelic variants
  3. They can cause merge errors and analysis issues

This tool should be run:
  - After downloading/converting data to PLINK format
  - Before phenotype generation
  - Before data partitioning/splitting

Methods:
  - PLINK2 (if available): Uses --max-alleles 2 flag (most reliable)
  - PLINK 1.9: Uses --freq to detect and --exclude to remove multiallelic variants

Example:
  python filter_multiallelic.py \\
      --plink-prefix data/genotypes \\
      --output-prefix data/genotypes_biallelic
        """
    )
    parser.add_argument(
        '--plink-prefix',
        type=str,
        required=True,
        help='Input PLINK prefix (path without .bed/.bim/.fam extension)'
    )
    parser.add_argument(
        '--output-prefix',
        type=str,
        required=True,
        help='Output PLINK prefix for filtered data'
    )
    parser.add_argument(
        '--keep-log',
        action='store_true',
        help='Keep PLINK log file (default: remove)'
    )
    
    args = parser.parse_args()
    
    # Validate input files exist
    plink_path = Path(args.plink_prefix)
    required_files = ['.bed', '.bim', '.fam']
    missing_files = []
    for ext in required_files:
        if not plink_path.with_suffix(ext).exists():
            missing_files.append(f"{args.plink_prefix}{ext}")
    
    if missing_files:
        logger.error(f"Missing required PLINK files: {missing_files}")
        sys.exit(1)
    
    # Create output directory if needed
    output_path = Path(args.output_prefix)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Filter multiallelic variants
    success = filter_multiallelic_variants(
        args.plink_prefix,
        args.output_prefix,
        keep_log=args.keep_log
    )
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
