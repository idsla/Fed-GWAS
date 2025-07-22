#!/usr/bin/env python3
# test_prg_masking.py

"""
Test script for PRG-MASKING secure aggregation implementation.
This script validates that the pairwise masking correctly aggregates data
while preserving privacy.
"""

import sys
import os
import numpy as np

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from server.prg_masking import PRGMaskingAggregator, ClientMaskingHelper

def test_prg_masking_correctness():
    """Test that PRG masking produces correct aggregation results."""
    print("=== Testing PRG-MASKING Correctness ===")
    
    num_clients = 3
    aggregator = PRGMaskingAggregator(num_clients)
    
    # Test data
    test_seeds = [12345, 67890, 54321]
    test_arrays = [
        np.array([1, 2, 3], dtype=np.float64),
        np.array([4, 5, 6], dtype=np.float64), 
        np.array([7, 8, 9], dtype=np.float64)
    ]
    
    # Expected results (plaintext aggregation)
    expected_seed_sum = sum(test_seeds)
    expected_array_sum = np.sum(test_arrays, axis=0)
    
    print(f"Original seeds: {test_seeds}")
    print(f"Expected seed sum: {expected_seed_sum}")
    print(f"Original arrays: {test_arrays}")
    print(f"Expected array sum: {expected_array_sum}")
    
    # Step 1: Simulate key exchange
    print("\n--- Step 1: Key Exchange ---")
    client_helpers = []
    dh_params = aggregator.get_dh_params()
    
    for client_id in range(num_clients):
        helper = ClientMaskingHelper(client_id, num_clients)
        public_key = helper.generate_dh_keypair(dh_params)
        aggregator.add_client_public_key(client_id, public_key)
        client_helpers.append(helper)
        print(f"Client {client_id} generated public key")
    
    assert aggregator.is_key_exchange_complete(), "Key exchange should be complete"
    all_public_keys = aggregator.get_all_public_keys()
    
    # Step 2: Clients compute shared secrets and mask data
    print("\n--- Step 2: Data Masking ---")
    masked_seeds = []
    masked_arrays = []
    
    for i, helper in enumerate(client_helpers):
        helper.compute_shared_secrets(all_public_keys, dh_params)
        
        # Mask seed
        masked_seed = helper.mask_data(test_seeds[i])
        masked_seeds.append(masked_seed)
        
        # Mask array
        masked_array = helper.mask_data(test_arrays[i])
        masked_arrays.append(masked_array)
        
        print(f"Client {i}: original seed {test_seeds[i]} -> masked {masked_seed:.2f}")
        print(f"Client {i}: original array {test_arrays[i]} -> masked {masked_array}")
    
    # Step 3: Server aggregation
    print("\n--- Step 3: Server Aggregation ---")
    aggregated_seed = sum(masked_seeds)
    aggregated_array = np.sum(masked_arrays, axis=0)
    
    print(f"Aggregated seed: {aggregated_seed:.2f}")
    print(f"Aggregated array: {aggregated_array}")
    
    # Step 4: Verify correctness
    print("\n--- Step 4: Verification ---")
    seed_error = abs(aggregated_seed - expected_seed_sum)
    array_error = np.abs(aggregated_array - expected_array_sum)
    
    print(f"Seed aggregation error: {seed_error:.2f}")
    print(f"Array aggregation error: {array_error}")
    
    # Allow small numerical errors due to floating point arithmetic
    assert seed_error < 1e-10, f"Seed aggregation error too large: {seed_error}"
    assert np.all(array_error < 1e-10), f"Array aggregation error too large: {array_error}"
    
    print("✅ PRG-MASKING correctness test PASSED!")
    return True

def test_privacy_properties():
    """Test that individual client data cannot be easily inferred."""
    print("\n=== Testing Privacy Properties ===")
    
    num_clients = 3
    aggregator = PRGMaskingAggregator(num_clients)
    
    # Test with sensitive data
    sensitive_seeds = [1000000, 2000000, 3000000]  # Large values
    
    # Key exchange
    client_helpers = []
    dh_params = aggregator.get_dh_params()
    
    for client_id in range(num_clients):
        helper = ClientMaskingHelper(client_id, num_clients)
        public_key = helper.generate_dh_keypair(dh_params)
        aggregator.add_client_public_key(client_id, public_key)
        client_helpers.append(helper)
    
    all_public_keys = aggregator.get_all_public_keys()
    
    # Mask data
    masked_seeds = []
    for i, helper in enumerate(client_helpers):
        helper.compute_shared_secrets(all_public_keys, dh_params)
        masked_seed = helper.mask_data(sensitive_seeds[i])
        masked_seeds.append(masked_seed)
        
        # Check that masked value is significantly different from original
        difference = abs(masked_seed - sensitive_seeds[i])
        print(f"Client {i}: original {sensitive_seeds[i]} -> masked {masked_seed:.2f} (diff: {difference:.2f})")
        assert difference > 1000, f"Masking difference too small for client {i}: {difference}"
    
    # Verify that individual values are well-hidden
    max_original = max(sensitive_seeds)
    masked_values_normalized = [abs(x) for x in masked_seeds]
    
    print(f"Max original value: {max_original}")
    print(f"Masked values range: {min(masked_values_normalized):.2f} - {max(masked_values_normalized):.2f}")
    
    # At least one masked value should be very different from any original
    min_distance_to_any_original = float('inf')
    for masked in masked_seeds:
        for original in sensitive_seeds:
            distance = abs(masked - original)
            min_distance_to_any_original = min(min_distance_to_any_original, distance)
    
    print(f"Minimum distance from any masked to any original: {min_distance_to_any_original:.2f}")
    assert min_distance_to_any_original > 100, "Some masked values too close to originals"
    
    print("✅ Privacy properties test PASSED!")
    return True

def test_multiple_rounds():
    """Test that masking works consistently across multiple rounds."""
    print("\n=== Testing Multiple Rounds ===")
    
    num_clients = 3
    num_rounds = 5
    
    for round_num in range(num_rounds):
        print(f"\n--- Round {round_num + 1} ---")
        
        aggregator = PRGMaskingAggregator(num_clients)
        
        # Different test data each round
        test_data = [round_num * 10 + i for i in range(num_clients)]
        expected_sum = sum(test_data)
        
        # Key exchange
        client_helpers = []
        dh_params = aggregator.get_dh_params()
        
        for client_id in range(num_clients):
            helper = ClientMaskingHelper(client_id, num_clients)
            public_key = helper.generate_dh_keypair(dh_params)
            aggregator.add_client_public_key(client_id, public_key)
            client_helpers.append(helper)
        
        all_public_keys = aggregator.get_all_public_keys()
        
        # Mask and aggregate
        masked_values = []
        for i, helper in enumerate(client_helpers):
            helper.compute_shared_secrets(all_public_keys, dh_params)
            masked_value = helper.mask_data(test_data[i])
            masked_values.append(masked_value)
        
        aggregated_value = sum(masked_values)
        error = abs(aggregated_value - expected_sum)
        
        print(f"Round {round_num + 1}: expected {expected_sum}, got {aggregated_value:.2f}, error {error:.2f}")
        assert error < 1e-10, f"Round {round_num + 1} failed with error {error}"
    
    print("✅ Multiple rounds test PASSED!")
    return True

def main():
    """Run all tests."""
    print("Starting PRG-MASKING Tests\n")
    
    try:
        test_prg_masking_correctness()
        test_privacy_properties()
        test_multiple_rounds()
        
        print("\n🎉 All PRG-MASKING tests PASSED!")
        print("\nPRG-MASKING implementation is ready for integration with Fed-GWAS pipeline.")
        
    except Exception as e:
        print(f"\n❌ Test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)