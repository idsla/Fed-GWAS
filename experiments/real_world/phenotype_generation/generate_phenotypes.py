#!/usr/bin/env python3
"""
Phenotype Generation Tool

Generate binary (case/control) phenotypes from PLINK genotype data.
Supports multiple models: additive, polygenic, and population-stratified.

Usage:
    python generate_phenotypes.py \
        --plink-prefix data/1000genomes \
        --output-fam data/1000genomes_pheno.fam \
        --model additive \
        --causal-fraction 0.01 \
        --seed 42
"""

import argparse
import logging
import json
from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed

from .models import AdditiveModel, PolygenicModel, PopulationStratifiedModel

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_genotypes(plink_prefix: str) -> np.ndarray:
    """
    Load genotypes from PLINK files.
    
    Args:
        plink_prefix: Path to PLINK files (without .bed/.bim/.fam extension)
        
    Returns:
        Genotype matrix (n_samples, n_snps) with values 0, 1, 2, or -9 (missing)
    """
    logger.info(f"Loading genotypes from {plink_prefix}")
    
    # Use pysnptools to read PLINK files
    try:
        snp_data = Bed(plink_prefix, count_A1=True)
        genotypes = snp_data.read().val
        
        # Convert to 0, 1, 2 format (pysnptools uses 0, 1, 2)
        # Missing values are already NaN, convert to -9
        genotypes = np.where(np.isnan(genotypes), -9, genotypes).astype(int)
        
        logger.info(f"Loaded {genotypes.shape[0]} samples and {genotypes.shape[1]} SNPs")
        return genotypes
    
    except Exception as e:
        logger.error(f"Failed to load genotypes: {e}")
        raise


def load_fam_file(plink_prefix: str) -> pd.DataFrame:
    """
    Load .fam file to get sample IDs.
    
    Args:
        plink_prefix: Path to PLINK files
        
    Returns:
        DataFrame with FID, IID, Father, Mother, Sex columns
    """
    fam_file = f"{plink_prefix}.fam"
    if not Path(fam_file).exists():
        raise FileNotFoundError(f"FAM file not found: {fam_file}")
    
    fam_df = pd.read_csv(
        fam_file,
        sep='\t',
        header=None,
        names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype']
    )
    return fam_df


def load_population_labels(fam_file: str, population_file: Optional[str] = None) -> Optional[np.ndarray]:
    """
    Load population labels if available.
    
    Args:
        fam_file: Path to .fam file
        population_file: Optional file with population labels (one per line, same order as .fam)
        
    Returns:
        Array of population labels or None
    """
    if population_file is None:
        return None
    
    if not Path(population_file).exists():
        logger.warning(f"Population file not found: {population_file}")
        return None
    
    labels = []
    with open(population_file, 'r') as f:
        for line in f:
            labels.append(line.strip())
    
    logger.info(f"Loaded {len(labels)} population labels")
    return np.array(labels)


def generate_phenotypes(
    plink_prefix: str,
    output_fam: str,
    model_type: str = 'additive',
    causal_fraction: float = 0.01,
    heritability: Optional[float] = None,
    effect_size_std: float = 0.1,
    population_effect_std: float = 0.2,
    noise_std: float = 1.0,
    intercept: float = 0.0,
    population_file: Optional[str] = None,
    seed: Optional[int] = None,
    metadata_file: Optional[str] = None
) -> None:
    """
    Generate phenotypes and update .fam file.
    
    Args:
        plink_prefix: Path to PLINK files (without extension)
        output_fam: Output .fam file path
        model_type: Model type ('additive', 'polygenic', 'population_stratified')
        causal_fraction: Fraction of causal SNPs
        heritability: Target heritability (for polygenic model)
        effect_size_std: Standard deviation of effect sizes
        population_effect_std: Standard deviation of population effects
        noise_std: Standard deviation of noise
        intercept: Intercept term
        population_file: Optional file with population labels
        seed: Random seed
        metadata_file: Optional file to save model metadata
    """
    # Load genotypes
    genotypes = load_genotypes(plink_prefix)
    
    # Load .fam file to get sample structure
    fam_df = load_fam_file(plink_prefix)
    
    # Load population labels if provided
    population_labels = load_population_labels(
        f"{plink_prefix}.fam",
        population_file
    )
    
    # Initialize model
    if model_type == 'additive':
        model = AdditiveModel(
            causal_fraction=causal_fraction,
            effect_size_std=effect_size_std,
            noise_std=noise_std,
            intercept=intercept,
            seed=seed
        )
        phenotypes, metadata = model.generate(genotypes)
    
    elif model_type == 'polygenic':
        if heritability is None:
            heritability = 0.3  # default
        model = PolygenicModel(
            causal_fraction=causal_fraction,
            heritability=heritability,
            intercept=intercept,
            seed=seed
        )
        phenotypes, metadata = model.generate(genotypes)
    
    elif model_type == 'population_stratified':
        model = PopulationStratifiedModel(
            causal_fraction=causal_fraction,
            effect_size_std=effect_size_std,
            population_effect_std=population_effect_std,
            noise_std=noise_std,
            intercept=intercept,
            seed=seed
        )
        phenotypes, metadata = model.generate(genotypes, population_labels)
    
    else:
        raise ValueError(f"Unknown model type: {model_type}")
    
    # Convert binary phenotypes to PLINK format (0->1 control, 1->2 case)
    plink_phenotypes = phenotypes + 1
    
    # Update .fam file with new phenotypes
    fam_df['Phenotype'] = plink_phenotypes
    
    # Save updated .fam file
    output_path = Path(output_fam)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fam_df.to_csv(output_path, sep='\t', index=False, header=False)
    
    logger.info(f"Generated phenotypes: {metadata['case_count']} cases, {metadata['control_count']} controls")
    logger.info(f"Case rate: {metadata['case_rate']:.3f}")
    logger.info(f"Saved updated .fam file to {output_fam}")
    
    # Save metadata if requested
    if metadata_file:
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        logger.info(f"Saved model metadata to {metadata_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate phenotypes from PLINK genotype data'
    )
    parser.add_argument(
        '--plink-prefix',
        type=str,
        required=True,
        help='Path to PLINK files (without .bed/.bim/.fam extension)'
    )
    parser.add_argument(
        '--output-fam',
        type=str,
        required=True,
        help='Output .fam file path'
    )
    parser.add_argument(
        '--model',
        type=str,
        choices=['additive', 'polygenic', 'population_stratified'],
        default='additive',
        help='Phenotype generation model (default: additive)'
    )
    parser.add_argument(
        '--causal-fraction',
        type=float,
        default=0.01,
        help='Fraction of SNPs that are causal (default: 0.01)'
    )
    parser.add_argument(
        '--heritability',
        type=float,
        default=None,
        help='Target heritability for polygenic model (default: 0.3)'
    )
    parser.add_argument(
        '--effect-size-std',
        type=float,
        default=0.1,
        help='Standard deviation of effect sizes (default: 0.1)'
    )
    parser.add_argument(
        '--population-effect-std',
        type=float,
        default=0.2,
        help='Standard deviation of population effects (default: 0.2)'
    )
    parser.add_argument(
        '--noise-std',
        type=float,
        default=1.0,
        help='Standard deviation of noise (default: 1.0)'
    )
    parser.add_argument(
        '--intercept',
        type=float,
        default=0.0,
        help='Intercept term (default: 0.0)'
    )
    parser.add_argument(
        '--population-file',
        type=str,
        default=None,
        help='Optional file with population labels (one per line)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=None,
        help='Random seed for reproducibility'
    )
    parser.add_argument(
        '--metadata-file',
        type=str,
        default=None,
        help='Optional file to save model metadata (JSON)'
    )
    
    args = parser.parse_args()
    
    generate_phenotypes(
        plink_prefix=args.plink_prefix,
        output_fam=args.output_fam,
        model_type=args.model,
        causal_fraction=args.causal_fraction,
        heritability=args.heritability,
        effect_size_std=args.effect_size_std,
        population_effect_std=args.population_effect_std,
        noise_std=args.noise_std,
        intercept=args.intercept,
        population_file=args.population_file,
        seed=args.seed,
        metadata_file=args.metadata_file
    )


if __name__ == '__main__':
    main()
