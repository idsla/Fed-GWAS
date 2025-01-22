# Kinship Analysis Module

## Overview

The **Kinship Analysis** module is designed to calculate kinship coefficients between individuals based on genetic data. It includes functionality to calculate the **KING kinship coefficient** and perform **incremental analysis** to estimate kinship relationships across multiple iterations. The module can handle large datasets and produce final kinship estimates that classify relationships into degrees (1st degree, 2nd degree, 3rd degree, or unrelated).



## Detailed User Guide

### 1. `calculate_king_coeff`

**Description**:  
This function calculates the KING kinship coefficient between two individuals using their genotype data.

**Parameters**:

- `genotype_i (np.ndarray)`: The genotype data for the first individual.
- `genotype_j (np.ndarray)`: The genotype data for the second individual.

**Returns**:

- `float`: The calculated KING kinship coefficient.

**Usage Example**:

```python
import numpy as np
from kinship_module import KinshipAnalyzer

kinship_analyzer = KinshipAnalyzer()

genotype_i = np.array([0, 1, 2, 0, 1])
genotype_j = np.array([0, 1, 2, 2, 1])

kinship_coefficient = kinship_analyzer.calculate_king_coeff(genotype_i, genotype_j)
print(f"KING kinship coefficient: {kinship_coefficient}")
```

---

### 2. `incremental_analysis`

**Description**:  
This function performs an incremental analysis to estimate kinship coefficients between individuals across multiple iterations. Each iteration uses an increasing number of SNPs to improve the estimation.

**Parameters**:

- `Da (np.ndarray)`: Genotype data for the first group of individuals.
- `Db (np.ndarray)`: Genotype data for the second group of individuals.
- `snps (np.ndarray)`: Indices of SNPs used in the analysis.
- `iterations (int)`: The number of iterations for incremental analysis.

**Returns**:

- `tuple`: A tuple containing the combined kinship coefficients (`combined_kinship`) and the kinship coefficients from each iteration (`kinship_history`).

**Usage Example**:

```python
import numpy as np
from kinship_module import KinshipAnalyzer

kinship_analyzer = KinshipAnalyzer()

Da = np.random.randint(0, 3, (10, 500))  # Simulated genotype data for 10 individuals
Db = np.random.randint(0, 3, (10, 500))  # Simulated genotype data for another 10 individuals
snps = np.arange(500)  # SNP indices
iterations = 10  # Number of iterations

combined_kinship, kinship_history = kinship_analyzer.incremental_analysis(Da, Db, snps, iterations)
print(f"Combined Kinship Coefficients: {combined_kinship}")
```

---

### 3. `final_kinship_estimate`

**Description**:  
This function provides the final classification of kinship relationships based on the kinship coefficient values. It classifies the relationships into 1st degree, 2nd degree, 3rd degree, or unrelated.

**Parameters**:

- `combined_kinship (dict)`: The final combined kinship coefficients after incremental analysis.

**Returns**:

- `dict`: A dictionary classifying the kinship relationships.

**Usage Example**:

```python
from kinship_module import KinshipAnalyzer

kinship_analyzer = KinshipAnalyzer()

# Assuming combined_kinship is obtained from incremental_analysis
final_estimates = kinship_analyzer.final_kinship_estimate(combined_kinship)
print(f"Final Kinship Estimates: {final_estimates}")
```

---

### 4. `KinshipLogger.monitor_kinship_history`

**Description**:  
This function logs the kinship history for each iteration of the incremental analysis, providing a detailed view of the changes in kinship estimates over time.

**Parameters**:

- `kinship_history (list)`: A list of dictionaries, each containing kinship coefficients for a particular iteration.

**Usage Example**:

```python
from kinship_module import KinshipLogger

# Assuming kinship_history is obtained from incremental_analysis
KinshipLogger.monitor_kinship_history(kinship_history)
```