"""Real-world data preparation adapters for packaged simulation examples."""

from __future__ import annotations

import random
import shutil
import subprocess
import urllib.request
from pathlib import Path

import typer


GENOMES_BASE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
GENOMES_VCF_NAME = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
GENOMES_VCF_INDEX_NAME = f"{GENOMES_VCF_NAME}.tbi"
GENOMES_PANEL_NAME = "integrated_call_samples_v3.20130502.ALL.panel"


def _download_if_missing(url: str, destination: Path) -> None:
    """Download a URL only when the destination file is absent.

    Args:
        url: Public HTTPS URL to retrieve.
        destination: Local file path inside the simulation project.

    Returns:
        None. Creates parent directories and writes the downloaded file.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        typer.echo(f"Already exists: {destination}")
        return
    typer.echo(f"Downloading {url}")
    urllib.request.urlretrieve(url, destination)


def _plink_binary() -> str:
    """Resolve a PLINK executable for real-data preparation.

    Returns:
        Executable name or absolute path.

    Raises:
        typer.BadParameter: If neither `plink` nor `plink2` is available.
    """
    for name in ("plink", "plink2"):
        found = shutil.which(name)
        if found:
            return found
    raise typer.BadParameter("PLINK is required to prepare 1000 Genomes data.")


def _run(cmd: list[str], project_dir: Path) -> None:
    """Run a preparation command in the project directory.

    Args:
        cmd: Command and arguments.
        project_dir: Simulation project root.

    Returns:
        None. Raises if the subprocess fails.
    """
    typer.echo("Running: " + " ".join(cmd))
    subprocess.run(cmd, cwd=project_dir, check=True)


def _write_binary_phenotypes(source_fam: Path, target_fam: Path, seed: int | None) -> None:
    """Copy a FAM file while replacing phenotypes with deterministic cases.

    Args:
        source_fam: Input PLINK `.fam` file.
        target_fam: Output `.fam` file with binary 1/2 phenotypes.
        seed: Optional seed for reproducible phenotype assignment.

    Returns:
        None. Writes `target_fam`.
    """
    rng = random.Random(seed)
    lines = []
    for line in source_fam.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) < 6:
            continue
        parts[5] = "2" if rng.random() < 0.5 else "1"
        lines.append("\t".join(parts))
    target_fam.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _sample_partitions(panel_path: Path) -> tuple[list[str], list[str]]:
    """Split 1000 Genomes sample ids into two population-based centers.

    Args:
        panel_path: Downloaded 1000 Genomes panel file.

    Returns:
        Tuple `(center_1_samples, center_2_samples)` with sample ids.
    """
    center_1: list[str] = []
    center_2: list[str] = []
    for line in panel_path.read_text(encoding="utf-8").splitlines():
        if not line or line.lower().startswith("sample"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        sample_id = parts[0]
        super_population = parts[2]
        if super_population == "EUR":
            center_1.append(sample_id)
        elif super_population in {"AFR", "AMR", "EAS", "SAS"}:
            center_2.append(sample_id)
    if not center_1 or not center_2:
        raise typer.BadParameter(f"Could not build two population partitions from {panel_path}.")
    return center_1, center_2


def _write_keep_file(path: Path, sample_ids: list[str]) -> None:
    """Write a PLINK keep file with FID/IID sample id pairs.

    Args:
        path: Destination keep file.
        sample_ids: Sample ids from the 1000 Genomes panel.

    Returns:
        None. Writes one `sample sample` row per id.
    """
    path.write_text("".join(f"{sample_id}\t{sample_id}\n" for sample_id in sample_ids), encoding="utf-8")


def prepare_1000genomes_chr22(project_dir: Path, seed: int | None = 42) -> None:
    """Download and prepare the packaged 1000 Genomes chr22 example.

    Args:
        project_dir: FedGWAS simulation project root.
        seed: Optional seed used for simple binary phenotype assignment.

    Returns:
        None. Writes PLINK triplets under `data/center_1` and `data/center_2`.

    Raises:
        typer.BadParameter: If PLINK is unavailable or sample partitions cannot
            be built from the downloaded panel file.
        subprocess.CalledProcessError: If a PLINK conversion command fails.
    """
    project_dir = project_dir.resolve()
    data_dir = project_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    vcf_path = data_dir / GENOMES_VCF_NAME
    panel_path = data_dir / GENOMES_PANEL_NAME

    _download_if_missing(f"{GENOMES_BASE_URL}/{GENOMES_VCF_NAME}", vcf_path)
    _download_if_missing(f"{GENOMES_BASE_URL}/{GENOMES_VCF_INDEX_NAME}", data_dir / GENOMES_VCF_INDEX_NAME)
    _download_if_missing(f"{GENOMES_BASE_URL}/{GENOMES_PANEL_NAME}", panel_path)

    plink = _plink_binary()
    raw_prefix = data_dir / "genotypes_raw"
    filtered_prefix = data_dir / "genotypes"
    pheno_prefix = data_dir / "genotypes_pheno"

    if not raw_prefix.with_suffix(".bed").exists():
        _run([plink, "--vcf", str(vcf_path), "--double-id", "--make-bed", "--out", str(raw_prefix)], project_dir)
    if not filtered_prefix.with_suffix(".bed").exists():
        _run(
            [
                plink,
                "--bfile",
                str(raw_prefix),
                "--biallelic-only",
                "strict",
                "--make-bed",
                "--out",
                str(filtered_prefix),
            ],
            project_dir,
        )

    for suffix in (".bed", ".bim"):
        shutil.copyfile(filtered_prefix.with_suffix(suffix), pheno_prefix.with_suffix(suffix))
    _write_binary_phenotypes(filtered_prefix.with_suffix(".fam"), pheno_prefix.with_suffix(".fam"), seed)

    center_1_samples, center_2_samples = _sample_partitions(panel_path)
    for center_id, samples in ((1, center_1_samples), (2, center_2_samples)):
        keep_file = data_dir / f"center_{center_id}_samples.txt"
        center_prefix = data_dir / f"center_{center_id}" / "genotypes"
        center_prefix.parent.mkdir(parents=True, exist_ok=True)
        _write_keep_file(keep_file, samples)
        _run(
            [
                plink,
                "--bfile",
                str(pheno_prefix),
                "--keep",
                str(keep_file),
                "--make-bed",
                "--out",
                str(center_prefix),
            ],
            project_dir,
        )
