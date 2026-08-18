#!/usr/bin/env python3
"""
Comprehensive test script for the visualization module.

This script tests all visualization functions with synthetic data to ensure
they work correctly in both interactive and headless (Docker) environments.
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import json

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from pipeline.visualization.plots_gwas import (
    make_gwas_plots_from_assoc_file,
    qq_plot,
    manhattan_plot,
    pvalue_hist,
)
from pipeline.visualization.plots_qc import (
    make_qc_plots,
    maf_from_genotype_counts,
    missing_rate_from_counts,
)
from pipeline.visualization.plots_king import (
    make_king_plots,
    parse_king_partial_results,
)
from pipeline.visualization.plots_system import (
    make_system_plots,
    plot_memory_usage_over_time,
    plot_cpu_usage_over_time,
    plot_network_usage,
    plot_resource_summary,
)
from pipeline.visualization.post_analysis import generate_post_analysis_report
from pipeline.visualization.hooks import (
    maybe_plot_client_global_qc,
    maybe_plot_client_local_lr,
    maybe_plot_client_king,
)


def create_test_assoc_file(output_dir: Path) -> Path:
    """Create a test PLINK .assoc.logistic file."""
    assoc_file = output_dir / "test.assoc.logistic"
    
    # Create realistic GWAS results
    n_snps = 1000
    chromosomes = np.random.choice(range(1, 23), n_snps)
    positions = np.random.randint(1, 100000000, n_snps)
    pvalues = np.random.beta(0.1, 10, n_snps)  # Most SNPs have low p-values
    
    # Add some significant hits
    significant_indices = np.random.choice(n_snps, 10, replace=False)
    pvalues[significant_indices] = np.random.uniform(1e-8, 1e-5, 10)
    
    data = {
        "CHR": chromosomes,
        "SNP": [f"rs{i}" for i in range(n_snps)],
        "BP": positions,
        "A1": np.random.choice(["A", "T", "G", "C"], n_snps),
        "A2": np.random.choice(["A", "T", "G", "C"], n_snps),
        "TEST": ["ADD"] * n_snps,
        "NMISS": np.random.randint(100, 500, n_snps),
        "OR": np.random.uniform(0.8, 1.2, n_snps),
        "SE": np.random.uniform(0.01, 0.1, n_snps),
        "L95": np.random.uniform(0.7, 1.1, n_snps),
        "U95": np.random.uniform(1.1, 1.3, n_snps),
        "STAT": np.random.normal(0, 1, n_snps),
        "P": pvalues,
    }
    
    df = pd.DataFrame(data)
    df.to_csv(assoc_file, sep="\t", index=False)
    return assoc_file


def create_test_qc_data():
    """Create test QC data (genotype counts and missingness)."""
    n_snps = 500
    counts = np.random.randint(0, 200, (n_snps, 3))
    missing = np.random.randint(0, 50, (n_snps, 2))
    return counts, missing


def create_test_king_data():
    """Create test KING results text."""
    n_pairs = 100
    lines = []
    for i in range(n_pairs):
        sample1 = f"sample_{i}"
        sample2 = f"sample_{(i+1) % n_pairs}"
        phi = np.random.uniform(-0.1, 0.5)
        n1 = np.random.randint(1000, 5000)
        lines.append(f"{sample1} {sample2} {phi} {n1}")
    return "\n".join(lines)


def create_test_metrics(output_dir: Path):
    """Create test performance metrics files."""
    # Stage metrics
    n_stages = 10
    n_clients = 3
    stages = ["key_exchange", "sync", "local_qc", "global_qc", "init_chunks",
              "iterative_king", "local_lr", "init_chunks_lr", "iterative_lr"]
    
    stage_data = []
    base_time = 1000000000.0
    
    for client_id in range(n_clients):
        current_time = base_time
        for stage_name in stages:
            duration = np.random.uniform(1, 30)
            stage_data.append({
                "stage_name": stage_name,
                "client_id": f"client_{client_id}",
                "start_time": current_time,
                "end_time": current_time + duration,
                "duration": duration,
                "memory_peak_mb": np.random.uniform(100, 1000),
                "cpu_percent": np.random.uniform(10, 80),
                "bytes_sent": np.random.randint(1000, 1000000),
                "bytes_received": np.random.randint(1000, 1000000),
                "success": True,
            })
            current_time += duration
    
    stage_df = pd.DataFrame(stage_data)
    stage_df.to_csv(output_dir / "stage_metrics.csv", index=False)
    
    # Client metrics
    client_data = []
    for client_id in range(n_clients):
        client_stages = stage_df[stage_df["client_id"] == f"client_{client_id}"]
        client_data.append({
            "client_id": f"client_{client_id}",
            "total_runtime": client_stages["duration"].sum(),
            "stages_completed": ",".join(client_stages["stage_name"].unique()),
            "stages_failed": "",
            "total_memory_peak_mb": client_stages["memory_peak_mb"].max(),
            "total_bytes_sent": client_stages["bytes_sent"].sum(),
            "total_bytes_received": client_stages["bytes_received"].sum(),
        })
    
    client_df = pd.DataFrame(client_data)
    client_df.to_csv(output_dir / "client_metrics.csv", index=False)
    
    # Network stats
    network_stats = []
    for i in range(50):
        network_stats.append({
            "timestamp": base_time + i * 10,
            "interface": "eth0",
            "bytes_sent": i * 10000,
            "bytes_recv": i * 15000,
            "packets_sent": i * 100,
            "packets_recv": i * 150,
            "errin": 0,
            "errout": 0,
            "dropin": 0,
            "dropout": 0,
        })
    
    with open(output_dir / "network_stats.json", "w") as f:
        json.dump(network_stats, f, indent=2)


def test_gwas_plots():
    """Test GWAS visualization functions."""
    print("Testing GWAS plots...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        assoc_file = create_test_assoc_file(tmpdir)
        
        # Test individual plots
        df = pd.read_csv(assoc_file, sep="\t")
        pvals = df["P"].values
        
        qq_path = qq_plot(pvals, tmpdir / "test_qq.png", title="Test QQ Plot")
        assert Path(qq_path).exists(), "QQ plot not created"
        print(f"  ✓ QQ plot: {qq_path}")
        
        manhattan_path = manhattan_plot(df, tmpdir / "test_manhattan.png", title="Test Manhattan")
        assert Path(manhattan_path).exists(), "Manhattan plot not created"
        print(f"  ✓ Manhattan plot: {manhattan_path}")
        
        hist_path = pvalue_hist(pvals, tmpdir / "test_hist.png", title="Test Histogram")
        assert Path(hist_path).exists(), "P-value histogram not created"
        print(f"  ✓ P-value histogram: {hist_path}")
        
        # Test combined function
        qq, manhattan, hist = make_gwas_plots_from_assoc_file(
            assoc_file, tmpdir, prefix="test_gwas"
        )
        assert all(p is not None for p in [qq, manhattan, hist]), "Some plots missing"
        print(f"  ✓ Combined GWAS plots function")
    print("  ✅ GWAS plots test passed\n")


def test_qc_plots():
    """Test QC visualization functions."""
    print("Testing QC plots...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        counts, missing = create_test_qc_data()
        
        maf_path, miss_path, hwe_path = make_qc_plots(
            counts, missing, tmpdir, prefix="test_qc", include_hwe=True
        )
        
        assert maf_path is not None and Path(maf_path).exists(), "MAF plot not created"
        assert miss_path is not None and Path(miss_path).exists(), "Missing plot not created"
        assert hwe_path is not None and Path(hwe_path).exists(), "HWE plot not created"
        
        print(f"  ✓ MAF histogram: {maf_path}")
        print(f"  ✓ Missing rate histogram: {miss_path}")
        print(f"  ✓ HWE p-value histogram: {hwe_path}")
    print("  ✅ QC plots test passed\n")


def test_king_plots():
    """Test KING visualization functions."""
    print("Testing KING plots...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        king_text = create_test_king_data()
        
        heatmap_path, dist_path, scatter_path = make_king_plots(
            king_text=king_text, out_dir=tmpdir, prefix="test_king"
        )
        
        assert heatmap_path is not None, "Heatmap not created"
        assert dist_path is not None, "Distribution plot not created"
        assert scatter_path is not None, "Scatter plot not created"
        
        print(f"  ✓ Kinship heatmap: {heatmap_path}")
        print(f"  ✓ Kinship distribution: {dist_path}")
        print(f"  ✓ Relationship scatter: {scatter_path}")
    print("  ✅ KING plots test passed\n")


def test_system_plots():
    """Test system statistics visualization."""
    print("Testing system plots...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        create_test_metrics(tmpdir)
        
        plots = make_system_plots(tmpdir, tmpdir / "system", title_prefix="Test")
        
        assert len(plots) > 0, "No system plots created"
        for name, path in plots.items():
            if path:
                assert Path(path).exists(), f"{name} plot not created"
                print(f"  ✓ {name}: {path}")
    print("  ✅ System plots test passed\n")


def test_client_hooks():
    """Test client-side visualization hooks."""
    print("Testing client hooks...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        log_dir = str(tmpdir)
        
        # Test QC hook
        counts, missing = create_test_qc_data()
        artifacts = maybe_plot_client_global_qc(
            client_id="test_client",
            log_dir=log_dir,
            counts_array=counts,
            missing_array=missing,
            viz_config={"enabled": True, "client_qc_plots": True},
        )
        assert artifacts is not None, "QC hook failed"
        print(f"  ✓ QC hook: {len(artifacts)} artifacts")
        
        # Test LR hook
        assoc_file = create_test_assoc_file(tmpdir)
        artifacts = maybe_plot_client_local_lr(
            client_id="test_client",
            log_dir=log_dir,
            assoc_file=str(assoc_file),
            viz_config={"enabled": True, "client_lr_plots": True},
        )
        assert artifacts is not None, "LR hook failed"
        print(f"  ✓ LR hook: {len(artifacts)} artifacts")
        
        # Test KING hook
        king_text = create_test_king_data()
        artifacts = maybe_plot_client_king(
            client_id="test_client",
            log_dir=log_dir,
            king_text=king_text,
            viz_config={"enabled": True, "client_king_plots": True},
        )
        assert artifacts is not None, "KING hook failed"
        print(f"  ✓ KING hook: {len(artifacts)} artifacts")
    print("  ✅ Client hooks test passed\n")


def test_post_analysis():
    """Test post-analysis report generation."""
    print("Testing post-analysis report...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create test data structure
        run_dir = tmpdir / "run"
        run_dir.mkdir()
        
        # Create server intermediate directory
        server_dir = run_dir / "server" / "intermediate" / "session123"
        server_dir.mkdir(parents=True)
        
        # Create test LR results
        assoc_file = create_test_assoc_file(server_dir)
        (server_dir / "lr_results.assoc.logistic").write_text(assoc_file.read_text())
        
        # Create metrics
        create_test_metrics(run_dir)
        
        # Generate report
        artifacts = generate_post_analysis_report(
            run_dir=run_dir,
            out_dir=run_dir / "viz",
            server_intermediate_dir=run_dir / "server" / "intermediate",
        )
        
        assert "report_md" in artifacts, "Report not generated"
        assert Path(artifacts["report_md"]).exists(), "Report file not created"
        
        print(f"  ✓ Report generated: {artifacts['report_md']}")
        print(f"  ✓ Total artifacts: {len(artifacts)}")
    print("  ✅ Post-analysis test passed\n")


def main():
    """Run all visualization tests."""
    print("=" * 60)
    print("Fed-GWAS Visualization Module Test Suite")
    print("=" * 60)
    print()
    
    # Set non-interactive backend for headless testing
    import matplotlib
    matplotlib.use("Agg")
    
    tests = [
        ("GWAS Plots", test_gwas_plots),
        ("QC Plots", test_qc_plots),
        ("KING Plots", test_king_plots),
        ("System Plots", test_system_plots),
        ("Client Hooks", test_client_hooks),
        ("Post-Analysis", test_post_analysis),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"  ❌ {test_name} test failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
            print()
    
    print("=" * 60)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("=" * 60)
    
    if failed == 0:
        print("✅ All visualization tests passed!")
        return 0
    else:
        print("❌ Some tests failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())

