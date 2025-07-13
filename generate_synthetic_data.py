#!/usr/bin/env python3
"""
Synthetic GWAS Data Generator

This script generates synthetic PLINK-format datasets at four scales:
- Tiny: 500 samples, 5,000 SNPs, 2 centers
- Small: 5,000 samples, 50,000 SNPs, 3 centers  
- Medium: 50,000 samples, 500,000 SNPs, 5 centers
- Large: 275,812 samples, 98,000,000 SNPs, 7 centers

Each scale is partitioned horizontally (by samples) across centers.
"""

import os
import numpy as np
import pandas as pd
import subprocess
import argparse
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SyntheticDataGenerator:
    def __init__(self, output_dir="simulated_data"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Dataset configurations
        self.scales = {
            'tiny': {
                'samples': 500,
                'snps': 5000,
                'centers': 2,
                'causal_fraction': 0.01,  # 1% causal SNPs
                'missing_rate': 0.02,     # 2% missing data
                'maf_range': (0.01, 0.5)
            },
            'small': {
                'samples': 5000,
                'snps': 50000,
                'centers': 3,
                'causal_fraction': 0.01,
                'missing_rate': 0.02,
                'maf_range': (0.01, 0.5)
            },
            'medium': {
                'samples': 50000,
                'snps': 500000,
                'centers': 5,
                'causal_fraction': 0.01,
                'missing_rate': 0.02,
                'maf_range': (0.01, 0.5)
            },
            'large': {
                'samples': 275812,
                'snps': 98000000,
                'centers': 7,
                'causal_fraction': 0.01,
                'missing_rate': 0.02,
                'maf_range': (0.01, 0.5)
            }
        }
    
    def generate_mafs(self, n_snps, maf_range=(0.01, 0.5)):
        """Generate minor allele frequencies for SNPs"""
        return np.random.uniform(maf_range[0], maf_range[1], n_snps)
    
    def generate_genotypes(self, n_samples, n_snps, mafs):
        """Generate genotypes based on MAFs"""
        genotypes = np.zeros((n_samples, n_snps), dtype=int)
        
        for i, maf in enumerate(mafs):
            # Sample genotypes from Binomial(2, MAF)
            genotypes[:, i] = np.random.binomial(2, maf, n_samples)
        
        return genotypes
    
    def introduce_missingness(self, genotypes, missing_rate=0.02):
        """Introduce random missing data"""
        n_samples, n_snps = genotypes.shape
        missing_mask = np.random.random((n_samples, n_snps)) < missing_rate
        
        # Use -9 for missing values (PLINK standard)
        genotypes_with_missing = genotypes.copy()
        genotypes_with_missing[missing_mask] = -9
        
        return genotypes_with_missing
    
    def generate_binary_phenotype(self, genotypes, causal_fraction=0.01):
        """Generate binary phenotype using logistic model"""
        n_samples, n_snps = genotypes.shape
        n_causal = int(n_snps * causal_fraction)
        
        # Randomly select causal SNPs
        causal_indices = np.random.choice(n_snps, n_causal, replace=False)
        
        # Generate effect sizes for causal SNPs
        effect_sizes = np.random.normal(0, 0.1, n_causal)
        
        # Calculate linear predictor
        linear_predictor = np.zeros(n_samples)
        for i, snp_idx in enumerate(causal_indices):
            # Handle missing values (treat as 0 for phenotype calculation)
            snp_data = genotypes[:, snp_idx]
            snp_data_clean = np.where(snp_data == -9, 0, snp_data)
            linear_predictor += effect_sizes[i] * snp_data_clean
        
        # Add noise
        linear_predictor += np.random.normal(0, 1, n_samples)
        
        # Convert to binary phenotype using logistic function
        probabilities = 1 / (1 + np.exp(-linear_predictor))
        phenotypes = np.random.binomial(1, probabilities)
        
        return phenotypes, causal_indices
    
    def create_plink_files(self, genotypes, phenotypes, mafs, scale_name):
        """Create PLINK-format files (.bed, .bim, .fam)"""
        n_samples, n_snps = genotypes.shape
        
        # Create output directory for this scale
        scale_dir = self.output_dir / scale_name
        scale_dir.mkdir(exist_ok=True)
        
        # Create .fam file (sample information)
        fam_data = []
        for i in range(n_samples):
            # FID, IID, Father, Mother, Sex, Phenotype
            fam_data.append([f"F{i+1:06d}", f"I{i+1:06d}", "0", "0", "1", phenotypes[i]])
        
        fam_df = pd.DataFrame(fam_data, columns=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])
        fam_file = scale_dir / f"{scale_name}.fam"
        fam_df.to_csv(fam_file, sep='\t', index=False, header=False)
        
        # Create .bim file (SNP information)
        bim_data = []
        for i in range(n_snps):
            # Chromosome, SNP ID, Genetic distance, Physical position, Allele1, Allele2
            bim_data.append([1, f"SNP{i+1:08d}", 0, i+1, "A", "T"])
        
        bim_df = pd.DataFrame(bim_data, columns=['Chr', 'SNP', 'CM', 'BP', 'A1', 'A2'])
        bim_file = scale_dir / f"{scale_name}.bim"
        bim_df.to_csv(bim_file, sep='\t', index=False, header=False)
        
        # Create .bed file (genotype data)
        bed_file = scale_dir / f"{scale_name}.bed"
        self._write_bed_file(genotypes, bed_file)
        
        logger.info(f"Created PLINK files for {scale_name}: {fam_file}, {bim_file}, {bed_file}")
        return scale_dir
    
    def _write_bed_file(self, genotypes, bed_file):
        """Write genotypes to PLINK .bed format"""
        n_samples, n_snps = genotypes.shape
        
        # PLINK .bed format: 2 bits per genotype
        # 00 = homozygous reference (0), 01 = missing (-9), 10 = heterozygous (1), 11 = homozygous alternate (2)
        
        with open(bed_file, 'wb') as f:
            # Write PLINK magic numbers
            f.write(b'\x6c\x1b\x01')
            
            # Write genotype data
            for snp_idx in range(n_snps):
                for sample_idx in range(0, n_samples, 4):
                    byte_val = 0
                    for bit_pos in range(4):
                        if sample_idx + bit_pos < n_samples:
                            genotype = genotypes[sample_idx + bit_pos, snp_idx]
                            if genotype == -9:  # missing
                                byte_val |= (1 << (bit_pos * 2))
                            elif genotype == 1:  # heterozygous
                                byte_val |= (2 << (bit_pos * 2))
                            elif genotype == 2:  # homozygous alternate
                                byte_val |= (3 << (bit_pos * 2))
                            # genotype == 0 (homozygous reference) stays as 0
                    f.write(bytes([byte_val]))
    
    def partition_samples(self, genotypes, phenotypes, n_centers, scale_name):
        """Partition samples across centers"""
        n_samples = genotypes.shape[0]
        
        # Create sample lists for each center
        for center in range(n_centers):
            center_dir = self.output_dir / scale_name / f"center_{center+1}"
            center_dir.mkdir(exist_ok=True)
            
            # Calculate sample indices for this center
            start_idx = (center * n_samples) // n_centers
            end_idx = ((center + 1) * n_samples) // n_centers
            sample_indices = list(range(start_idx, end_idx))
            
            # Create sample list file
            sample_list_file = center_dir / "samples.txt"
            with open(sample_list_file, 'w') as f:
                for idx in sample_indices:
                    f.write(f"F{idx+1:06d}\tI{idx+1:06d}\n")
            
            # Use PLINK to create center-specific files
            input_prefix = self.output_dir / scale_name / scale_name
            output_prefix = center_dir / f"{scale_name}_center_{center+1}"
            
            plink_cmd = [
                "plink",
                "--bfile", str(input_prefix),
                "--keep", str(sample_list_file),
                "--make-bed",
                "--out", str(output_prefix)
            ]
            
            try:
                result = subprocess.run(plink_cmd, capture_output=True, text=True, check=True)
                logger.info(f"Created center {center+1} files for {scale_name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"PLINK command failed for center {center+1}: {e}")
                logger.error(f"STDOUT: {e.stdout}")
                logger.error(f"STDERR: {e.stderr}")
    
    def generate_scale(self, scale_name):
        """Generate data for a specific scale"""
        logger.info(f"Generating {scale_name} scale dataset...")
        
        config = self.scales[scale_name]
        n_samples = config['samples']
        n_snps = config['snps']
        n_centers = config['centers']
        
        # Generate MAFs
        logger.info(f"Generating MAFs for {n_snps} SNPs...")
        mafs = self.generate_mafs(n_snps, config['maf_range'])
        
        # Generate genotypes
        logger.info(f"Generating genotypes for {n_samples} samples...")
        genotypes = self.generate_genotypes(n_samples, n_snps, mafs)
        
        # Introduce missingness
        logger.info("Introducing missing data...")
        genotypes = self.introduce_missingness(genotypes, config['missing_rate'])
        
        # Generate binary phenotype
        logger.info("Generating binary phenotype...")
        phenotypes, causal_snps = self.generate_binary_phenotype(genotypes, config['causal_fraction'])
        
        # Create PLINK files
        logger.info("Creating PLINK files...")
        scale_dir = self.create_plink_files(genotypes, phenotypes, mafs, scale_name)
        
        # Partition samples across centers
        logger.info(f"Partitioning samples across {n_centers} centers...")
        self.partition_samples(genotypes, phenotypes, n_centers, scale_name)
        
        # Save metadata
        metadata = {
            'scale': scale_name,
            'n_samples': n_samples,
            'n_snps': n_snps,
            'n_centers': n_centers,
            'causal_snps': causal_snps.tolist(),
            'case_count': int(np.sum(phenotypes)),
            'control_count': int(n_samples - np.sum(phenotypes))
        }
        
        metadata_file = scale_dir / "metadata.json"
        import json
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Completed {scale_name} scale generation!")
        return scale_dir
    
    def generate_all_scales(self):
        """Generate data for all scales"""
        for scale_name in self.scales.keys():
            try:
                self.generate_scale(scale_name)
            except Exception as e:
                logger.error(f"Failed to generate {scale_name} scale: {e}")
                continue

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic GWAS datasets")
    parser.add_argument("--scale", choices=['tiny', 'small', 'medium', 'large', 'all'], 
                       default='all', help="Scale to generate")
    parser.add_argument("--output-dir", default="simulated_data", 
                       help="Output directory")
    
    args = parser.parse_args()
    
    generator = SyntheticDataGenerator(args.output_dir)
    
    if args.scale == 'all':
        generator.generate_all_scales()
    else:
        generator.generate_scale(args.scale)
    
    logger.info("Synthetic data generation completed!")

if __name__ == "__main__":
    main() 