#!/usr/bin/env python3
"""
Test runner for synthetic data generation

Run this script to execute all tests for the synthetic data generation pipeline.
"""

import subprocess
import sys
import os

def main():
    """Run the test suite"""
    # Get the directory containing this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    test_dir = os.path.join(script_dir, "tests")
    
    print("Running synthetic data generation tests...")
    print(f"Test directory: {test_dir}")
    print("=" * 60)
    
    # Run pytest
    cmd = [
        sys.executable, "-m", "pytest",
        os.path.join(test_dir, "test_synthetic_data.py"),
        "-v",
        "--tb=short"
    ]
    
    try:
        result = subprocess.run(cmd, check=True)
        print("\n" + "=" * 60)
        print("✅ All tests passed!")
        return 0
    except subprocess.CalledProcessError as e:
        print("\n" + "=" * 60)
        print("❌ Some tests failed!")
        return e.returncode

if __name__ == "__main__":
    sys.exit(main()) 