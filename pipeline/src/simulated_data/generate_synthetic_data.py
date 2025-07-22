#!/usr/bin/env python3
"""
Synthetic GWAS Data Generator

This script generates synthetic PLINK-format datasets at four scales:
- Tiny: 500 samples, 5,000 SNPs, 2 centers
- Small: 5,000 samples, 50,000 SNPs, 3 centers  
- Medium: 50,000 samples, 500,000 SNPs, 5 centers
- Large: 275,812 samples, 98,000,000 SNPs, 7 centers

Each scale is partitioned horizontally (by samples) across centers.
Includes family relationships with first-, second-, and third-degree relatives.
"""

import os
import numpy as np
import pandas as pd
import subprocess
import argparse
from pathlib import Path
import logging
import random
import json
import csv
import shutil

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
                'maf_range': (0.01, 0.5),
                'relative_fraction': 0.1,  # 15% of samples are relatives
                'relative_distribution': {'first': 0.4, 'second': 0.35, 'third': 0.25}  # Distribution of relative types
            },
            'small': {
                'samples': 5000,
                'snps': 50000,
                'centers': 3,
                'causal_fraction': 0.01,
                'missing_rate': 0.02,
                'maf_range': (0.01, 0.5),
                'relative_fraction': 0.1,
                'relative_distribution': {'first': 0.4, 'second': 0.35, 'third': 0.25}
            },
            'medium': {
                'samples': 50000,
                'snps': 500000,
                'centers': 5,
                'causal_fraction': 0.01,
                'missing_rate': 0.02,
                'maf_range': (0.01, 0.5),
                'relative_fraction': 0.1,
                'relative_distribution': {'first': 0.4, 'second': 0.35, 'third': 0.25}
            },
            'large': {
                'samples': 275812,
                'snps': 98000000,
                'centers': 7,
                'causal_fraction': 0.01,
                'missing_rate': 0.02,
                'maf_range': (0.01, 0.5),
                'relative_fraction': 0.1,
                'relative_distribution': {'first': 0.4, 'second': 0.35, 'third': 0.25}
            }
        }
    
    def generate_mafs(self, n_snps, maf_range=(0.01, 0.5)):
        """Generate minor allele frequencies for SNPs"""
        return np.random.uniform(maf_range[0], maf_range[1], n_snps)
    
    def generate_family_structure(self, n_samples, relative_fraction=0.1, relative_distribution=None):
        """
        Generate family structure with relatives
        
        Args:
            n_samples: Total number of samples
            relative_fraction: Fraction of samples that are relatives
            relative_distribution: Dict with fractions for first, second, third degree relatives
        
        Returns:
            family_info: Dict with family assignments and relationships
        """
        if relative_distribution is None:
            relative_distribution = {'first': 0.5, 'second': 0.35, 'third': 0.15}
        
        n_relatives = int(n_samples * relative_fraction)
        n_unrelated = n_samples - n_relatives
        
        # Start with unrelated individuals
        family_info = {
            'family_ids': np.zeros(n_samples, dtype=int),
            'relationship_types': np.full(n_samples, 'unrelated', dtype=object),
            'parent_ids': np.full(n_samples, -1, dtype=int),
            'spouse_ids': np.full(n_samples, -1, dtype=int)
        }
        
        # Assign family IDs to unrelated individuals
        family_info['family_ids'][:n_unrelated] = np.arange(n_unrelated)
        
        # Generate relatives
        relative_idx = n_unrelated
        family_id = n_unrelated
        
        # Calculate number of each type of relative
        n_first = int(n_relatives * relative_distribution['first'])
        n_second = int(n_relatives * relative_distribution['second'])
        n_third = int(n_relatives * relative_distribution['third'])
        
        # Generate first-degree relatives (parent-child, siblings)
        for i in range(n_first):
            if relative_idx >= n_samples:
                break
                
            # Randomly select a base individual from unrelated samples
            base_idx = np.random.randint(0, n_unrelated)
            
            # 50% chance of parent-child, 50% chance of sibling
            if np.random.random() < 0.5:
                # Parent-child relationship
                family_info['family_ids'][relative_idx] = family_info['family_ids'][base_idx]
                family_info['relationship_types'][relative_idx] = 'parent_child'
                family_info['parent_ids'][relative_idx] = base_idx
            else:
                # Sibling relationship
                family_info['family_ids'][relative_idx] = family_info['family_ids'][base_idx]
                family_info['relationship_types'][relative_idx] = 'sibling'
                family_info['parent_ids'][relative_idx] = family_info['parent_ids'][base_idx]
            
            relative_idx += 1
        
        # Generate second-degree relatives (grandparent-grandchild, uncle/aunt-niece/nephew)
        for i in range(n_second):
            if relative_idx >= n_samples:
                break
                
            base_idx = np.random.randint(0, n_unrelated)
            family_info['family_ids'][relative_idx] = family_info['family_ids'][base_idx]
            family_info['relationship_types'][relative_idx] = 'second_degree'
            relative_idx += 1
        
        # Generate third-degree relatives (cousins, etc.)
        for i in range(n_third):
            if relative_idx >= n_samples:
                break
                
            base_idx = np.random.randint(0, n_unrelated)
            family_info['family_ids'][relative_idx] = family_info['family_ids'][base_idx]
            family_info['relationship_types'][relative_idx] = 'third_degree'
            relative_idx += 1
        
        return family_info
    
    def generate_genotypes(self, n_samples, n_snps, mafs, family_info=None):
        """
        Generate genotypes based on MAFs, accounting for family relationships
        
        Args:
            n_samples: Number of samples
            n_snps: Number of SNPs
            mafs: Minor allele frequencies
            family_info: Family structure information (optional)
        
        Returns:
            genotypes: Genotype matrix
        """
        genotypes = np.zeros((n_samples, n_snps), dtype=int)
        
        # Generate base genotypes for all individuals
        for i, maf in enumerate(mafs):
            # Sample genotypes from Binomial(2, MAF) for unrelated individuals
            genotypes[:, i] = np.random.binomial(2, maf, n_samples)
        
        # If family structure is provided, modify genotypes for relatives
        if family_info is not None:
            genotypes = self._adjust_genotypes_for_relatedness(genotypes, family_info, mafs)
        
        return genotypes
    
    def _adjust_genotypes_for_relatedness(self, genotypes, family_info, mafs):
        """
        Adjust genotypes to reflect family relationships
        
        Args:
            genotypes: Base genotype matrix
            family_info: Family structure information
            mafs: Minor allele frequencies
        
        Returns:
            adjusted_genotypes: Genotype matrix with relatedness
        """
        n_samples, n_snps = genotypes.shape
        adjusted_genotypes = genotypes.copy()
        
        # Relationship coefficients (proportion of shared alleles)
        relationship_coeffs = {
            'parent_child': 0.5,      # Parent-child: 50% shared alleles
            'sibling': 0.5,           # Siblings: 50% shared alleles
            'second_degree': 0.25,    # Second-degree: 25% shared alleles
            'third_degree': 0.125     # Third-degree: 12.5% shared alleles
        }
        
        for i in range(n_samples):
            relationship = family_info['relationship_types'][i]
            
            if relationship == 'unrelated':
                continue
            
            # Find the base individual this person is related to
            if family_info['parent_ids'][i] >= 0:
                base_idx = family_info['parent_ids'][i]
            else:
                # For other relationships, randomly select a base individual from same family
                family_id = family_info['family_ids'][i]
                family_members = np.where(family_info['family_ids'] == family_id)[0]
                base_candidates = [j for j in family_members if j != i and family_info['relationship_types'][j] == 'unrelated']
                if base_candidates:
                    base_idx = np.random.choice(base_candidates)
                else:
                    continue
            
            # Adjust genotypes based on relationship
            coeff = relationship_coeffs.get(relationship, 0.0)
            
            for snp_idx in range(n_snps):
                base_genotype = genotypes[base_idx, snp_idx]
                
                if base_genotype == -9:  # Missing in base individual
                    continue
                
                # With probability coeff, inherit allele from base individual
                if np.random.random() < coeff:
                    # Inherit one allele from base individual
                    if base_genotype == 0:  # Homozygous reference
                        inherited_alleles = 0
                    elif base_genotype == 1:  # Heterozygous
                        inherited_alleles = np.random.choice([0, 1])
                    else:  # Homozygous alternate
                        inherited_alleles = 1
                    
                    # Generate second allele based on MAF
                    maf = mafs[snp_idx]
                    second_allele = np.random.binomial(1, maf)
                    
                    # Combine alleles
                    adjusted_genotypes[i, snp_idx] = inherited_alleles + second_allele
        
        return adjusted_genotypes
    
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
    
    def create_plink_files(self, genotypes, phenotypes, mafs, scale_name, family_info=None):
        """Create PLINK-format files (.bed, .bim, .fam)"""
        n_samples, n_snps = genotypes.shape
        
        # Create output directory for this scale
        scale_dir = self.output_dir / scale_name
        scale_dir.mkdir(exist_ok=True)
        
        # Create .fam file (sample information)
        fam_data = []
        for i in range(n_samples):
            # FID, IID, Father, Mother, Sex, Phenotype
            if family_info is not None and family_info['parent_ids'][i] >= 0:
                # This is a relative with a parent
                parent_id = family_info['parent_ids'][i]
                father_id = f"F{parent_id+1:06d}" if family_info['relationship_types'][i] == 'parent_child' else "0"
                mother_id = f"F{parent_id+1:06d}" if family_info['relationship_types'][i] == 'parent_child' else "0"
            else:
                father_id = "0"
                mother_id = "0"
            
            fam_data.append([f"F{i+1:06d}", f"I{i+1:06d}", father_id, mother_id, "1", phenotypes[i]])
        
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
    
    def partition_samples(self, genotypes, phenotypes, n_centers, scale_name, partition_strategy="even"):
        """
        Partition samples across centers using different strategies
        
        Args:
            genotypes: Genotype matrix
            phenotypes: Phenotype array
            n_centers: Number of centers
            scale_name: Scale name for output
            partition_strategy: "even", "uneven", or "skewed"
        """
        n_samples = genotypes.shape[0]
        
        if partition_strategy == "even":
            # Even Sample Split - equal samples per center, balanced phenotype
            sample_indices_per_center = self._even_partition(n_samples, n_centers, phenotypes)
        elif partition_strategy == "uneven":
            # Uneven Sample Sizes - different sample counts, balanced phenotype
            sample_indices_per_center = self._uneven_partition(n_samples, n_centers, phenotypes)
        elif partition_strategy == "skewed":
            # Uneven Size + Phenotype Skew - different sample counts and phenotype ratios
            sample_indices_per_center = self._skewed_partition(n_samples, n_centers, phenotypes)
        else:
            raise ValueError(f"Unknown partition strategy: {partition_strategy}")
        
        # Create center-specific files
        all_assigned = set()
        for center in range(n_centers):
            center_dir = self.output_dir / scale_name / f"center_{center+1}"
            center_dir.mkdir(exist_ok=True)
            
            sample_indices = sample_indices_per_center[center]
            
            # Create sample list file
            sample_list_file = center_dir / "samples.txt"
            with open(sample_list_file, 'w') as f:
                for idx in sample_indices:
                    fid = f"F{idx+1:06d}"
                    iid = f"I{idx+1:06d}"
                    f.write(f"{fid}\t{iid}\n")
                    all_assigned.add((fid, iid))
            
            # --- Per-center sample ID mapping ---
            id_map_file = center_dir / f'id_map_center_{center+1}.csv'
            with open(id_map_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['original_index', 'FID', 'IID'])
                for idx in sample_indices:
                    writer.writerow([idx, f'F{idx+1:06d}', f'I{idx+1:06d}'])
            # --- End per-center sample ID mapping ---

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
                center_phenotypes = phenotypes[sample_indices]
                case_ratio = float(np.mean(center_phenotypes))
                logger.info(f"Created center {center+1} files for {scale_name}: {len(sample_indices)} samples, {case_ratio:.1%} cases")
                # --- Per-center metadata generation ---
                # Load relatives_info.csv if exists
                rel_info_path = self.output_dir / scale_name / "relatives_info.csv"
                relatives_in_center = []
                if rel_info_path.exists():
                    rel_df = pd.read_csv(rel_info_path)
                    # Build set of sample FIDs in this center
                    sample_fids = set([f"I{idx+1:06d}" for idx in sample_indices])
                    # Filter relatives where both SampleID and AncestryID are in this center
                    rels = rel_df[(rel_df['SampleID'].isin(sample_fids)) & (rel_df['AncestryID'].isin(sample_fids))]
                    relatives_in_center = rels.to_dict(orient='records')
                # Build metadata dict
                metadata = {
                    "center": center+1,
                    "sample_count": int(len(sample_indices)),
                    "case_count": int(np.sum(center_phenotypes)),
                    "control_count": int(len(sample_indices) - np.sum(center_phenotypes)),
                    "case_ratio": case_ratio,
                    "sample_ids": [f"F{idx+1:06d}" for idx in sample_indices],
                    "relatives": relatives_in_center
                }
                # Write metadata to file
                metadata_file = center_dir / f"metadata_center_{center+1}.json"
                with open(metadata_file, 'w') as mf:
                    json.dump(metadata, mf, indent=2)
                # --- End per-center metadata ---
            except subprocess.CalledProcessError as e:
                logger.error(f"PLINK command failed for center {center+1}: {e}")
                logger.error(f"STDOUT: {e.stdout}")
                logger.error(f"STDERR: {e.stderr}")
        # --- Write assigned_samples.txt in main dataset directory ---
        assigned_samples_file = self.output_dir / scale_name / "assigned_samples.txt"
        with open(assigned_samples_file, 'w') as f:
            for fid, iid in sorted(all_assigned):
                f.write(f"{fid}\t{iid}\n")
        # --- End assigned_samples.txt ---
    
    def _even_partition(self, n_samples, n_centers, phenotypes):
        """Even sample split with balanced phenotype distribution. Ensures all samples are assigned."""
        samples_per_center = n_samples // n_centers
        remainder = n_samples % n_centers
        
        # Separate cases and controls
        case_indices = np.where(phenotypes == 1)[0].tolist()
        control_indices = np.where(phenotypes == 0)[0].tolist()
        
        sample_indices_per_center = []
        for center in range(n_centers):
            center_samples = samples_per_center + (1 if center < remainder else 0)
            cases_needed = center_samples // 2
            controls_needed = center_samples - cases_needed
            center_cases = []
            center_controls = []
            # Assign as many cases as possible
            for _ in range(cases_needed):
                if case_indices:
                    center_cases.append(case_indices.pop())
                elif control_indices:
                    center_controls.append(control_indices.pop())
            # Assign as many controls as possible
            for _ in range(controls_needed):
                if control_indices:
                    center_controls.append(control_indices.pop())
                elif case_indices:
                    center_cases.append(case_indices.pop())
            center_indices = np.array(center_cases + center_controls)
            np.random.shuffle(center_indices)
            sample_indices_per_center.append(center_indices)
        # If any samples remain, assign them round-robin
        leftovers = case_indices + control_indices
        for i, idx in enumerate(leftovers):
            sample_indices_per_center[i % n_centers] = np.append(sample_indices_per_center[i % n_centers], idx)
        return sample_indices_per_center

    def _uneven_partition(self, n_samples, n_centers, phenotypes):
        """Hybrid: Assign as close as possible to ratios, then distribute leftovers round-robin."""
        if n_centers == 2:
            ratios = [0.6, 0.4]
        elif n_centers == 3:
            ratios = [0.5, 0.3, 0.2]
        elif n_centers == 5:
            ratios = [0.4, 0.25, 0.2, 0.1, 0.05]
        elif n_centers == 7:
            ratios = [0.3, 0.2, 0.15, 0.15, 0.1, 0.05, 0.05]
        else:
            ratios = [1.0 / n_centers] * n_centers
        case_indices = np.where(phenotypes == 1)[0].tolist()
        control_indices = np.where(phenotypes == 0)[0].tolist()
        sample_indices_per_center = []
        total_assigned = 0
        leftovers = []
        for center in range(n_centers):
            center_samples = int(n_samples * ratios[center])
            if center == n_centers - 1:
                center_samples = n_samples - total_assigned  # assign all remaining
            total_assigned += center_samples
            cases_needed = center_samples // 2
            controls_needed = center_samples - cases_needed
            center_cases = []
            center_controls = []
            for _ in range(cases_needed):
                if case_indices:
                    center_cases.append(case_indices.pop())
                elif control_indices:
                    leftovers.append(None)
            for _ in range(controls_needed):
                if control_indices:
                    center_controls.append(control_indices.pop())
                elif case_indices:
                    leftovers.append(None)
            center_indices = np.array(center_cases + center_controls)
            np.random.shuffle(center_indices)
            sample_indices_per_center.append(center_indices)
        leftovers = case_indices + control_indices
        for i, idx in enumerate(leftovers):
            sample_indices_per_center[i % n_centers] = np.append(sample_indices_per_center[i % n_centers], idx)
        # Log final counts and ratios
        for center, indices in enumerate(sample_indices_per_center):
            center_phenotypes = phenotypes[indices]
            case_ratio = np.mean(center_phenotypes)
            print(f"[Hybrid Partition] Center {center+1}: {len(indices)} samples, case ratio: {case_ratio:.2f}")
        return sample_indices_per_center

    def _skewed_partition(self, n_samples, n_centers, phenotypes):
        """Hybrid: Assign as close as possible to ratios, then distribute leftovers round-robin."""
        if n_centers == 3:
            size_ratios = [0.6, 0.3, 0.1]
            case_ratios = [0.5, 0.3, 0.1]
        elif n_centers == 5:
            size_ratios = [0.4, 0.25, 0.2, 0.1, 0.05]
            case_ratios = [0.5, 0.4, 0.3, 0.2, 0.1]
        elif n_centers == 7:
            size_ratios = [0.3, 0.2, 0.15, 0.15, 0.1, 0.05, 0.05]
            case_ratios = [0.5, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15]
        else:
            size_ratios = [1.0 / n_centers] * n_centers
            case_ratios = [0.5] * n_centers
        case_indices = np.where(phenotypes == 1)[0].tolist()
        control_indices = np.where(phenotypes == 0)[0].tolist()
        sample_indices_per_center = []
        total_assigned = 0
        leftovers = []
        for center in range(n_centers):
            center_samples = int(n_samples * size_ratios[center])
            if center == n_centers - 1:
                center_samples = n_samples - total_assigned
            total_assigned += center_samples
            case_ratio = case_ratios[center]
            cases_needed = int(center_samples * case_ratio)
            controls_needed = center_samples - cases_needed
            center_cases = []
            center_controls = []
            for _ in range(cases_needed):
                if case_indices:
                    center_cases.append(case_indices.pop())
                elif control_indices:
                    leftovers.append(None)
            for _ in range(controls_needed):
                if control_indices:
                    center_controls.append(control_indices.pop())
                elif case_indices:
                    leftovers.append(None)
            center_indices = np.array(center_cases + center_controls)
            np.random.shuffle(center_indices)
            sample_indices_per_center.append(center_indices)
        leftovers = case_indices + control_indices
        for i, idx in enumerate(leftovers):
            sample_indices_per_center[i % n_centers] = np.append(sample_indices_per_center[i % n_centers], idx)
        # Log final counts and ratios
        for center, indices in enumerate(sample_indices_per_center):
            center_phenotypes = phenotypes[indices]
            case_ratio = np.mean(center_phenotypes)
            print(f"[Hybrid Partition] Center {center+1}: {len(indices)} samples, case ratio: {case_ratio:.2f}")
        return sample_indices_per_center
    
    def generate_scale(self, scale_name, partition_strategy="even", seed=None):
        """Generate data for a specific scale using Mendelian relatives logic. Optionally fix random seed."""
        logger.info(f"Generating {scale_name} scale dataset with {partition_strategy} partitioning...")
        
        config = self.scales[scale_name]
        n_samples = config['samples']
        n_snps = config['snps']
        n_centers = config['centers']
        n_relatives = int(n_samples * config['relative_fraction'])
        n_unrelated = n_samples - n_relatives
        
        # Optionally set the random seed for reproducibility
        if seed is not None:
            import random
            import numpy as np
            random.seed(seed)
            np.random.seed(seed)
        
        # Generate MAFs
        logger.info(f"Generating MAFs for {n_snps} SNPs...")
        mafs = self.generate_mafs(n_snps, config['maf_range'])
        
        # Generate unrelated individuals
        logger.info(f"Generating {n_unrelated} unrelated individuals...")
        unrelated_genotypes = np.zeros((n_snps, n_unrelated), dtype=int)
        for i, maf in enumerate(mafs):
            unrelated_genotypes[i, :] = np.random.binomial(2, maf, n_unrelated)
        
        # Generate relatives using Mendelian inheritance
        logger.info(f"Generating {n_relatives} relatives using Mendelian inheritance...")
        Gdata_full, rel_ids_df = self.create_relatives(unrelated_genotypes, n=n_relatives)
        # Gdata_full: DataFrame (n_snps, n_unrelated + n_relatives)
        all_sample_ids = list(Gdata_full.columns)
        genotypes = Gdata_full.to_numpy().T  # shape: (n_samples, n_snps)
        
        # Generate binary phenotype
        logger.info("Generating binary phenotype...")
        phenotypes, causal_snps = self.generate_binary_phenotype(genotypes, config['causal_fraction'])
        
        # Create family_info for .fam and metadata
        family_info = {
            'family_ids': np.arange(n_samples),
            'relationship_types': np.array(['unrelated'] * n_unrelated + list(rel_ids_df['Relatedness'].values)),
            'parent_ids': np.array([-1] * n_unrelated + [-1] * n_relatives),  # Not tracked in this logic
            'spouse_ids': np.array([-1] * n_samples)
        }
        
        # Create PLINK files with new genotypes and family info
        logger.info("Creating PLINK files...")
        scale_dir = self.create_plink_files(genotypes, phenotypes, mafs, scale_name, family_info)
        
        # Partition samples across centers
        logger.info(f"Partitioning samples across {n_centers} centers using {partition_strategy} strategy...")
        self.partition_samples(genotypes, phenotypes, n_centers, scale_name, partition_strategy)
        
        # Save metadata including family information
        metadata = {
            'scale': scale_name,
            'partition_strategy': partition_strategy,
            'n_samples': n_samples,
            'n_snps': n_snps,
            'n_centers': n_centers,
            'causal_snps': causal_snps.tolist(),
            'case_count': int(np.sum(phenotypes)),
            'control_count': int(n_samples - np.sum(phenotypes)),
            'family_info': {
                'relative_fraction': config['relative_fraction'],
                'n_relatives': n_relatives,
                'n_unrelated': n_unrelated
            },
            'random_seed': seed
        }
        
        metadata_file = scale_dir / "metadata.json"
        import json
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # Save detailed relatives info
        rel_ids_file = scale_dir / "relatives_info.csv"
        rel_ids_df.to_csv(rel_ids_file, index=False)
        
        logger.info(f"Completed {scale_name} scale generation!")
        logger.info(f"Family structure: {n_relatives} relatives, {n_unrelated} unrelated")
        return scale_dir
    
    def generate_all_scales(self):
        """Generate data for all scales"""
        for scale_name in self.scales.keys():
            try:
                self.generate_scale(scale_name)
            except Exception as e:
                logger.error(f"Failed to generate {scale_name} scale: {e}")
                continue

    def add_child(self, df, childID, fatherID, motherID):
        '''
        Add child SNPs data to the dataframe using Mendelian inheritance
        df: (n_snps, n_users) - pandas DataFrame
        childID: str - child ID
        fatherID: str - father ID
        motherID: str - mother ID
        '''
        fatherCol = df[fatherID].to_numpy()
        motherCol = df[motherID].to_numpy()
        childCol = np.zeros_like(fatherCol)
        for i in range(len(fatherCol)):
            p = random.uniform(0, 1)
            if fatherCol[i] == 0 and motherCol[i] == 0:
                childCol[i] = 0
            elif (fatherCol[i] == 0 and motherCol[i] == 1) or (fatherCol[i] == 1 and motherCol[i] == 0):
                childCol[i] = 0 if p > 0.5 else 1
            elif (fatherCol[i] == 0 and motherCol[i] == 2) or (fatherCol[i] == 2 and motherCol[i] == 0):
                childCol[i] = 1
            elif fatherCol[i] == 1 and motherCol[i] == 1:
                if p < 0.333:
                    childCol[i] = 0
                elif p < 0.666:
                    childCol[i] = 1
                else:
                    childCol[i] = 2
            elif (fatherCol[i] == 1 and motherCol[i] == 2) or (fatherCol[i] == 2 and motherCol[i] == 1):
                childCol[i] = 1 if p > 0.5 else 2
            elif fatherCol[i] == 2 and motherCol[i] == 2:
                childCol[i] = 2
        df[childID] = childCol
        return df

    def create_relatives(self, Gdata, n=40):
        '''
        Create synthetic relatives using Mendelian inheritance
        Gdata: (n_snps, n_users) - numpy array
        n: int - number of relatives
        '''
        rel_ids = []
        Gdata_n = pd.DataFrame(Gdata)
        Gdata_n.columns = [f"I{i+1:06d}" for i in range(Gdata.shape[1])]
        parent_ids = np.random.choice(Gdata_n.columns, size=n, replace=False)
        relative_mapping = {}
        # Step 1: 1st generation
        F_id = parent_ids[:n//2]
        M_id = parent_ids[n//2:]
        child_ids_1C = [f'1C_{i}' for i in range(n//2)]
        for c_id, f_id, m_id in zip(child_ids_1C, F_id, M_id):
            Gdata_n = self.add_child(Gdata_n, c_id, f_id, m_id)
            relative_mapping[c_id] = [f_id, m_id]
            for ances_id in relative_mapping[c_id]:
                rel_ids.append([c_id, ances_id, 'first-degree'])
        # Step 2: 2nd generation
        np.random.shuffle(child_ids_1C)
        group_1CF = child_ids_1C[:n//4]
        group_1CM = child_ids_1C[n//4:]
        child_ids_2C = [f'2C_{i}' for i in range(n//4)]
        for c_id, f_id, m_id in zip(child_ids_2C, group_1CF, group_1CM):
            Gdata_n = self.add_child(Gdata_n, c_id, f_id, m_id)
            f_ances_ids = relative_mapping[f_id]
            m_ances_ids = relative_mapping[m_id]
            relative_mapping[c_id] = f_ances_ids + m_ances_ids
            for ansces_id in relative_mapping[c_id]:
                rel_ids.append([c_id, ansces_id, 'second-degree'])
        # Step 3: 3rd generation
        np.random.shuffle(child_ids_2C)
        group_2CF = child_ids_2C[:n//8]
        group_2CM = child_ids_2C[n//8:]
        child_ids_3C = [f'3C_{i}' for i in range(n//8)]
        for c_id, f_id, m_id in zip(child_ids_3C, group_2CF, group_2CM):
            Gdata_n = self.add_child(Gdata_n, c_id, f_id, m_id)
            f_ances_ids = relative_mapping[f_id]
            m_ances_ids = relative_mapping[m_id]
            for ansces_id in f_ances_ids + m_ances_ids:
                rel_ids.append([c_id, ansces_id, 'third-degree'])
        rel_ids_df = pd.DataFrame(rel_ids, columns=["SampleID", "AncestryID", "Relatedness"])
        return Gdata_n, rel_ids_df

    def calculate_king_kinship(self, genotype_i, genotype_j):
        n11 = np.sum((genotype_i == 1) & (genotype_j == 1))
        n02 = np.sum((genotype_i == 2) & (genotype_j == 0))
        n20 = np.sum((genotype_i == 0) & (genotype_j == 2))
        n1_s = np.sum(genotype_i == 1)
        s_1 = np.sum(genotype_j == 1)
        if n1_s == 0:
            return 0
        phi_ij = (2 * n11 - 4 * (n02 + n20) - n1_s + s_1) / (4 * n1_s)
        if phi_ij > 0.5:
            return 0.5
        if phi_ij < 0:
            return 0
        return phi_ij

    def generate_ground_truths(self, scale_name):
        """Run full QC+association pipeline for ground truth"""
        scale_dir = self.output_dir / scale_name
        centralized_dir = scale_dir / "centralized"
        centralized_dir.mkdir(exist_ok=True)
        bfile = scale_dir / scale_name
        plink = "plink"
        # Step 1: Initial QC - missingness, HWE
        qc1_prefix = centralized_dir / "qc1"
        qc1_cmd = [
            plink, "--bfile", str(bfile),
            "--mind", "0.05", "--geno", "0.05", "--hwe", "1e-6",
            "--make-bed", "--out", str(qc1_prefix)
        ]
        try:
            subprocess.run(qc1_cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"PLINK QC1 failed: {e}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
            return
        # Step 2: Kinship/relatedness on QCed data
        genome_cmd = [plink, "--bfile", str(qc1_prefix), "--genome", "--out", str(qc1_prefix)]
        try:
            subprocess.run(genome_cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"PLINK genome failed: {e}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
            return
        # Step 3: Remove relateds (PI_HAT > 0.2)
        genome_file = f"{qc1_prefix}.genome"
        fam_file = f"{qc1_prefix}.fam"
        unrelated_file = centralized_dir / "unrelated.txt"
        try:
            fam_df = pd.read_csv(fam_file, sep='\s+', header=None)
            fam_ids = set(fam_df[1].tolist())
            if pd.io.common.file_exists(genome_file):
                genome_df = pd.read_csv(genome_file, delim_whitespace=True)
                related = set()
                for _, row in genome_df.iterrows():
                    if row['PI_HAT'] > 0.2:
                        related.add(row['IID1'])
                        related.add(row['IID2'])
                unrelated = fam_ids - related
            else:
                unrelated = fam_ids
            unrelated_fam = fam_df[fam_df[1].isin(unrelated)]
            with open(unrelated_file, 'w') as f:
                for _, row in unrelated_fam.iterrows():
                    f.write(f"{row[0]} {row[1]}\n")
        except Exception as e:
            logger.error(f"Error processing relateds: {e}")
            return
        # Step 4: Create final unrelated, QCed dataset
        qc2_prefix = centralized_dir / "qc2"
        qc2_cmd = [
            plink, "--bfile", str(qc1_prefix), "--keep", str(unrelated_file),
            "--make-bed", "--out", str(qc2_prefix)
        ]
        try:
            subprocess.run(qc2_cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"PLINK QC2 (unrelated) failed: {e}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
            return
        # Step 5: Final association analysis
        assoc_cmd = [plink, "--bfile", str(qc2_prefix), "--assoc", "--out", str(centralized_dir / "assoc_final")]
        try:
            subprocess.run(assoc_cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"PLINK assoc_final failed: {e}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
            return
        logger.info("Centralized QC+association ground truth pipeline completed.")

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic GWAS datasets")
    parser.add_argument("--scale", choices=['tiny', 'small', 'medium', 'large', 'all'], 
                       default='all', help="Scale to generate")
    parser.add_argument("--partition-strategy", choices=['even', 'uneven', 'skewed'], 
                       default='even', help="Partition strategy: even, uneven, or skewed")
    parser.add_argument("--output-dir", default="simulated_data", 
                       help="Output directory")
    parser.add_argument("--ground-truths", action="store_true", help="Generate ground truth stats/results after generation")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    args = parser.parse_args()
    if args.seed is not None:
        import random
        import numpy as np
        random.seed(args.seed)
        np.random.seed(args.seed)
    generator = SyntheticDataGenerator(args.output_dir)
    
    if args.scale == 'all':
        # Generate all scales with the specified partition strategy
        for scale_name in generator.scales.keys():
            try:
                generator.generate_scale(scale_name, args.partition_strategy, seed=args.seed)
                if args.ground_truths:
                    generator.generate_ground_truths(scale_name)
            except Exception as e:
                logger.error(f"Failed to generate {scale_name} scale: {e}")
                continue
    else:
        generator.generate_scale(args.scale, args.partition_strategy, seed=args.seed)
        if args.ground_truths:
            generator.generate_ground_truths(args.scale)
    
    logger.info("Synthetic data generation completed!")

if __name__ == "__main__":
    main()
