import numpy as np

class QCUtils:

    @staticmethod
    def snphwe(obs_hets: int, obs_hom1: int, obs_hom2: int) -> float:
        """
        Perform Hardy-Weinberg Equilibrium (HWE) test using exact p-value calculation.

        Parameters:
        obs_hets (int): Number of heterozygous genotypes.
        obs_hom1 (int): Number of homozygous genotypes for the first allele.
        obs_hom2 (int): Number of homozygous genotypes for the second allele.

        Returns:
        float: The p-value for the Hardy-Weinberg equilibrium test.
        """
        if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
            raise ValueError("Observed counts must be non-negative integers.")
        
        # Total number of genotypes
        N = obs_hom1 + obs_hom2 + obs_hets
        
        # Identify rare and common homozygotes
        obs_homr = min(obs_hom1, obs_hom2)
        obs_homc = max(obs_hom1, obs_hom2)
        
        # Number of rare allele copies
        rare = obs_homr * 2 + obs_hets
        
        # Initialize probability array
        probs = np.zeros(rare + 1)
        
        # Find the midpoint of the distribution
        mid = rare * (2 * N - rare) // (2 * N)
        if (mid % 2) != (rare % 2):
            mid += 1
        
        probs[mid] = 1.0
        mysum = 1.0
        
        # Calculate probabilities from midpoint downwards
        curr_hets = mid
        curr_homr = (rare - mid) // 2
        curr_homc = N - curr_hets - curr_homr
        
        while curr_hets >= 2:
            probs[curr_hets - 2] = probs[curr_hets] * curr_hets * (curr_hets - 1) / (4 * (curr_homr + 1) * (curr_homc + 1))
            mysum += probs[curr_hets - 2]
            curr_hets -= 2
            curr_homr += 1
            curr_homc += 1
        
        # Calculate probabilities from midpoint upwards
        curr_hets = mid
        curr_homr = (rare - mid) // 2
        curr_homc = N - curr_hets - curr_homr
        
        while curr_hets <= rare - 2:
            probs[curr_hets + 2] = probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
            mysum += probs[curr_hets + 2]
            curr_hets += 2
            curr_homr -= 1
            curr_homc -= 1
        
        # P-value calculation
        target = probs[obs_hets]
        p_value = min(1.0, np.sum(probs[probs <= target]) / mysum)
        return p_value