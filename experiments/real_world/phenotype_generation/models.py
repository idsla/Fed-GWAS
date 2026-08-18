#!/usr/bin/env python3
"""
Phenotype Generation Models

This module provides models for generating phenotypes from genotypes.
Used for creating realistic case/control phenotypes on real genotype data
(e.g., 1000 Genomes Project) for federated GWAS evaluation.
"""

import numpy as np
from typing import List, Tuple, Optional, Dict
from scipy.special import expit  # logistic function


class PhenotypeModel:
    """Base class for phenotype generation models"""
    
    def __init__(self, seed: Optional[int] = None):
        """
        Initialize the model.
        
        Args:
            seed: Random seed for reproducibility
        """
        if seed is not None:
            np.random.seed(seed)
        self.seed = seed
    
    def generate(self, genotypes: np.ndarray, **kwargs) -> Tuple[np.ndarray, Dict]:
        """
        Generate phenotypes from genotypes.
        
        Args:
            genotypes: Genotype matrix (n_samples, n_snps) with values 0, 1, 2, or -9 (missing)
            **kwargs: Model-specific parameters
            
        Returns:
            Tuple of (phenotypes, metadata)
            - phenotypes: Binary phenotype array (0=control, 1=case)
            - metadata: Dictionary with model parameters and causal SNPs
        """
        raise NotImplementedError


class AdditiveModel(PhenotypeModel):
    """
    Additive logistic model for binary phenotype generation.
    
    Model: P(Y=1) = logistic(β₀ + Σᵢ βᵢGᵢ + ε)
    
    Where:
    - Gᵢ is the genotype at SNP i (0, 1, or 2)
    - βᵢ is the effect size for SNP i
    - ε is random noise
    """
    
    def __init__(
        self,
        causal_fraction: float = 0.01,
        effect_size_mean: float = 0.0,
        effect_size_std: float = 0.1,
        noise_std: float = 1.0,
        intercept: float = 0.0,
        seed: Optional[int] = None
    ):
        """
        Initialize additive model.
        
        Args:
            causal_fraction: Fraction of SNPs that are causal (default: 0.01 = 1%)
            effect_size_mean: Mean of effect size distribution (default: 0.0)
            effect_size_std: Standard deviation of effect size distribution (default: 0.1)
            noise_std: Standard deviation of random noise (default: 1.0)
            intercept: Intercept term β₀ (default: 0.0)
            seed: Random seed for reproducibility
        """
        super().__init__(seed)
        self.causal_fraction = causal_fraction
        self.effect_size_mean = effect_size_mean
        self.effect_size_std = effect_size_std
        self.noise_std = noise_std
        self.intercept = intercept
    
    def generate(
        self,
        genotypes: np.ndarray,
        causal_snps: Optional[List[int]] = None
    ) -> Tuple[np.ndarray, Dict]:
        """
        Generate binary phenotypes using additive model.
        
        Args:
            genotypes: Genotype matrix (n_samples, n_snps)
            causal_snps: Optional list of causal SNP indices. If None, randomly selected.
            
        Returns:
            Tuple of (phenotypes, metadata)
        """
        n_samples, n_snps = genotypes.shape
        
        # Select causal SNPs if not provided
        if causal_snps is None:
            n_causal = int(n_snps * self.causal_fraction)
            causal_snps = np.random.choice(n_snps, n_causal, replace=False).tolist()
        
        # Generate effect sizes for causal SNPs
        effect_sizes = np.random.normal(
            self.effect_size_mean,
            self.effect_size_std,
            len(causal_snps)
        )
        
        # Calculate linear predictor
        linear_predictor = np.full(n_samples, self.intercept)
        
        for i, snp_idx in enumerate(causal_snps):
            # Handle missing values (treat as 0 for phenotype calculation)
            snp_data = genotypes[:, snp_idx]
            snp_data_clean = np.where(snp_data == -9, 0, snp_data)
            linear_predictor += effect_sizes[i] * snp_data_clean
        
        # Add noise
        linear_predictor += np.random.normal(0, self.noise_std, n_samples)
        
        # Convert to binary phenotype using logistic function
        probabilities = expit(linear_predictor)
        phenotypes = np.random.binomial(1, probabilities)
        
        metadata = {
            'model': 'additive',
            'causal_snps': causal_snps,
            'effect_sizes': effect_sizes.tolist(),
            'causal_fraction': self.causal_fraction,
            'effect_size_mean': self.effect_size_mean,
            'effect_size_std': self.effect_size_std,
            'noise_std': self.noise_std,
            'intercept': self.intercept,
            'seed': self.seed,
            'case_count': int(np.sum(phenotypes)),
            'control_count': int(n_samples - np.sum(phenotypes)),
            'case_rate': float(np.mean(phenotypes))
        }
        
        return phenotypes, metadata


class PolygenicModel(PhenotypeModel):
    """
    Polygenic model with multiple causal variants.
    
    Similar to additive model but with explicit control over heritability.
    """
    
    def __init__(
        self,
        causal_fraction: float = 0.01,
        heritability: float = 0.3,
        intercept: float = 0.0,
        seed: Optional[int] = None
    ):
        """
        Initialize polygenic model.
        
        Args:
            causal_fraction: Fraction of SNPs that are causal
            heritability: Target heritability (h²) - proportion of variance explained by genetics
            intercept: Intercept term
            seed: Random seed for reproducibility
        """
        super().__init__(seed)
        self.causal_fraction = causal_fraction
        self.heritability = heritability
        self.intercept = intercept
    
    def generate(
        self,
        genotypes: np.ndarray,
        causal_snps: Optional[List[int]] = None
    ) -> Tuple[np.ndarray, Dict]:
        """
        Generate binary phenotypes using polygenic model.
        
        Args:
            genotypes: Genotype matrix (n_samples, n_snps)
            causal_snps: Optional list of causal SNP indices
            
        Returns:
            Tuple of (phenotypes, metadata)
        """
        n_samples, n_snps = genotypes.shape
        
        # Select causal SNPs if not provided
        if causal_snps is None:
            n_causal = int(n_snps * self.causal_fraction)
            causal_snps = np.random.choice(n_snps, n_causal, replace=False).tolist()
        
        # Calculate genetic component
        genetic_component = np.zeros(n_samples)
        effect_sizes = []
        
        for snp_idx in causal_snps:
            snp_data = genotypes[:, snp_idx]
            snp_data_clean = np.where(snp_data == -9, 0, snp_data)
            
            # Standardize genotypes
            snp_mean = np.mean(snp_data_clean)
            snp_std = np.std(snp_data_clean)
            if snp_std > 0:
                snp_data_std = (snp_data_clean - snp_mean) / snp_std
            else:
                snp_data_std = snp_data_clean
            
            # Generate effect size
            effect_size = np.random.normal(0, 1)
            effect_sizes.append(effect_size)
            genetic_component += effect_size * snp_data_std
        
        # Scale genetic component to achieve target heritability
        genetic_var = np.var(genetic_component)
        if genetic_var > 0:
            # Calculate required environmental variance
            # h² = var(G) / (var(G) + var(E))
            # var(E) = var(G) * (1 - h²) / h²
            env_var = genetic_var * (1 - self.heritability) / self.heritability
            env_std = np.sqrt(env_var)
        else:
            env_std = 1.0
        
        # Add environmental component
        environmental_component = np.random.normal(0, env_std, n_samples)
        linear_predictor = self.intercept + genetic_component + environmental_component
        
        # Convert to binary phenotype
        probabilities = expit(linear_predictor)
        phenotypes = np.random.binomial(1, probabilities)
        
        metadata = {
            'model': 'polygenic',
            'causal_snps': causal_snps,
            'effect_sizes': effect_sizes,
            'causal_fraction': self.causal_fraction,
            'heritability': self.heritability,
            'intercept': self.intercept,
            'seed': self.seed,
            'case_count': int(np.sum(phenotypes)),
            'control_count': int(n_samples - np.sum(phenotypes)),
            'case_rate': float(np.mean(phenotypes))
        }
        
        return phenotypes, metadata


class PopulationStratifiedModel(PhenotypeModel):
    """
    Model with population stratification effects.
    
    Adds population-specific effects to simulate population structure.
    Useful for testing federated GWAS on multi-population datasets.
    """
    
    def __init__(
        self,
        causal_fraction: float = 0.01,
        effect_size_std: float = 0.1,
        population_effect_std: float = 0.2,
        noise_std: float = 1.0,
        intercept: float = 0.0,
        seed: Optional[int] = None
    ):
        """
        Initialize population-stratified model.
        
        Args:
            causal_fraction: Fraction of SNPs that are causal
            effect_size_std: Standard deviation of SNP effect sizes
            population_effect_std: Standard deviation of population-specific effects
            noise_std: Standard deviation of random noise
            intercept: Intercept term
            seed: Random seed for reproducibility
        """
        super().__init__(seed)
        self.causal_fraction = causal_fraction
        self.effect_size_std = effect_size_std
        self.population_effect_std = population_effect_std
        self.noise_std = noise_std
        self.intercept = intercept
    
    def generate(
        self,
        genotypes: np.ndarray,
        population_labels: Optional[np.ndarray] = None,
        causal_snps: Optional[List[int]] = None
    ) -> Tuple[np.ndarray, Dict]:
        """
        Generate binary phenotypes with population stratification.
        
        Args:
            genotypes: Genotype matrix (n_samples, n_snps)
            population_labels: Array of population labels (e.g., 'EUR', 'AFR', 'ASN')
            causal_snps: Optional list of causal SNP indices
            
        Returns:
            Tuple of (phenotypes, metadata)
        """
        n_samples, n_snps = genotypes.shape
        
        # Select causal SNPs if not provided
        if causal_snps is None:
            n_causal = int(n_snps * self.causal_fraction)
            causal_snps = np.random.choice(n_snps, n_causal, replace=False).tolist()
        
        # Generate effect sizes
        effect_sizes = np.random.normal(0, self.effect_size_std, len(causal_snps))
        
        # Calculate genetic component
        linear_predictor = np.full(n_samples, self.intercept)
        
        for i, snp_idx in enumerate(causal_snps):
            snp_data = genotypes[:, snp_idx]
            snp_data_clean = np.where(snp_data == -9, 0, snp_data)
            linear_predictor += effect_sizes[i] * snp_data_clean
        
        # Add population-specific effects if labels provided
        population_effects = {}
        if population_labels is not None:
            unique_populations = np.unique(population_labels)
            for pop in unique_populations:
                pop_effect = np.random.normal(0, self.population_effect_std)
                population_effects[pop] = float(pop_effect)
                pop_mask = population_labels == pop
                linear_predictor[pop_mask] += pop_effect
        else:
            # If no labels, treat all as one population
            unique_populations = ['ALL']
            population_effects['ALL'] = 0.0
        
        # Add noise
        linear_predictor += np.random.normal(0, self.noise_std, n_samples)
        
        # Convert to binary phenotype
        probabilities = expit(linear_predictor)
        phenotypes = np.random.binomial(1, probabilities)
        
        metadata = {
            'model': 'population_stratified',
            'causal_snps': causal_snps,
            'effect_sizes': effect_sizes.tolist(),
            'population_effects': population_effects,
            'causal_fraction': self.causal_fraction,
            'effect_size_std': self.effect_size_std,
            'population_effect_std': self.population_effect_std,
            'noise_std': self.noise_std,
            'intercept': self.intercept,
            'seed': self.seed,
            'case_count': int(np.sum(phenotypes)),
            'control_count': int(n_samples - np.sum(phenotypes)),
            'case_rate': float(np.mean(phenotypes))
        }
        
        return phenotypes, metadata
