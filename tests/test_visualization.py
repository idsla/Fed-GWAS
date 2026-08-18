from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.visualization.plots_gwas import make_gwas_plots_from_assoc_file
from pipeline.visualization.plots_qc import make_qc_plots


def test_make_qc_plots(tmp_path: Path) -> None:
    # counts: [N_AA, N_Aa, N_aa]
    counts = np.array(
        [
            [10, 5, 0],
            [2, 8, 10],
            [20, 0, 0],
        ],
        dtype=np.int64,
    )
    # missing: [N_obs, N_miss]
    missing = np.array(
        [
            [100, 0],
            [90, 10],
            [50, 50],
        ],
        dtype=np.int64,
    )

    maf_path, miss_path = make_qc_plots(counts, missing, tmp_path, prefix="qc_test", title_prefix="test")
    # Either plot may be None if data is empty; here both should be present
    assert maf_path is not None and Path(maf_path).exists()
    assert miss_path is not None and Path(miss_path).exists()


def test_make_gwas_plots_from_assoc_file(tmp_path: Path) -> None:
    # Minimal assoc.logistic-like data
    df = pd.DataFrame(
        {
            "CHR": [1, 1, 2, 2],
            "SNP": ["rs1", "rs2", "rs3", "rs4"],
            "BP": [100, 200, 300, 400],
            "P": [0.5, 1e-4, 0.2, 0.9],
        }
    )
    assoc = tmp_path / "lr_results.assoc.logistic"
    df.to_csv(assoc, sep="\t", index=False)

    qq, manhattan, hist = make_gwas_plots_from_assoc_file(assoc, tmp_path, prefix="gwas_test", title_prefix="test")
    assert qq is not None and Path(qq).exists()
    assert hist is not None and Path(hist).exists()
    # manhattan can be None if parsing fails; should succeed here
    assert manhattan is not None and Path(manhattan).exists()


