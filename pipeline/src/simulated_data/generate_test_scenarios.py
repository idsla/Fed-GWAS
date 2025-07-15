#!/usr/bin/env python3
"""
Generate Test Scenarios Script

This script generates all the test scenarios outlined in the specification:
- A1, A2: Even sample split, balanced phenotype
- B1, B2, B3: Small scale with different partitioning strategies
- C1, C2, C3: Medium scale with different partitioning strategies
- D: Large scale with SF-GWAS partition strategy
"""

import subprocess
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_generation(scale, partition_strategy, output_suffix=""):
    """Run the synthetic data generation for a specific scenario"""
    output_dir = f"simulated_data{output_suffix}"
    
    cmd = [
        "python", "pipeline/src/simulated_data/generate_synthetic_data.py",
        "--scale", scale,
        "--partition-strategy", partition_strategy,
        "--output-dir", output_dir
    ]
    
    logger.info(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Successfully generated {scale} with {partition_strategy} partitioning")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to generate {scale} with {partition_strategy} partitioning")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        return False

def main():
    """Generate all test scenarios"""
    
    # Define scenarios based on the specification
    scenarios = [
        # A1: 1kGP - Even Sample Split, balanced phenotype (5 clients)
        # Note: This would require real 1kGP data, so we'll skip for now
        
        # A2: Sim-Tiny - Even Sample Split, balanced phenotype (2 clients)
        ("tiny", "even", "_A2"),
        
        # B1: Sim-Small - Even Sample Split, balanced phenotype (3 clients)
        ("small", "even", "_B1"),
        
        # B2: Sim-Small - Uneven Sample Sizes, balanced phenotype (3 clients)
        ("small", "uneven", "_B2"),
        
        # B3: Sim-Small - Uneven Size + Phenotype Skew (3 clients)
        ("small", "skewed", "_B3"),
        
        # C1: Sim-Medium - Even Sample Split, balanced phenotype (5 clients)
        ("medium", "even", "_C1"),
        
        # C2: Sim-Medium - Uneven Sample Sizes, balanced phenotype (5 clients)
        ("medium", "uneven", "_C2"),
        
        # C3: Sim-Medium - Uneven Size + Phenotype Skew (5 clients)
        ("medium", "skewed", "_C3"),
        
        # D: Sim-Large - SF-GWAS partition strategy (7 clients)
        # Using skewed strategy for SF-GWAS-like partitioning
        ("large", "skewed", "_D"),
    ]
    
    logger.info("Starting generation of all test scenarios...")
    
    successful_scenarios = []
    failed_scenarios = []
    
    for scale, partition_strategy, output_suffix in scenarios:
        logger.info(f"\n{'='*60}")
        logger.info(f"Generating Scenario: {scale} with {partition_strategy} partitioning")
        logger.info(f"{'='*60}")
        
        success = run_generation(scale, partition_strategy, output_suffix)
        
        if success:
            successful_scenarios.append((scale, partition_strategy, output_suffix))
        else:
            failed_scenarios.append((scale, partition_strategy, output_suffix))
    
    # Summary
    logger.info(f"\n{'='*60}")
    logger.info("GENERATION SUMMARY")
    logger.info(f"{'='*60}")
    logger.info(f"Successful scenarios: {len(successful_scenarios)}")
    logger.info(f"Failed scenarios: {len(failed_scenarios)}")
    
    if successful_scenarios:
        logger.info("\nSuccessful scenarios:")
        for scale, partition_strategy, output_suffix in successful_scenarios:
            logger.info(f"  - {scale} ({partition_strategy}) -> simulated_data{output_suffix}")
    
    if failed_scenarios:
        logger.info("\nFailed scenarios:")
        for scale, partition_strategy, output_suffix in failed_scenarios:
            logger.info(f"  - {scale} ({partition_strategy}) -> simulated_data{output_suffix}")
    
    logger.info("\nGeneration completed!")

if __name__ == "__main__":
    main() 