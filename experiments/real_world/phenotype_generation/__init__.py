"""
Phenotype Generation Tools

Tools for generating phenotypes from real genotype data (e.g., 1000 Genomes Project).
"""

from .models import (
    PhenotypeModel,
    AdditiveModel,
    PolygenicModel,
    PopulationStratifiedModel
)

__all__ = [
    'PhenotypeModel',
    'AdditiveModel',
    'PolygenicModel',
    'PopulationStratifiedModel',
]
