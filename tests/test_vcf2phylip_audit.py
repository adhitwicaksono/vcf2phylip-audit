import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "vcf2phylip.py"
DATA = ROOT / "tests" / "data"


def load_module():
    spec = importlib.util.spec_from_file_location("vcf2phylip", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_fasta(path):
    records = {}
    current = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            current = line[1:]
            records[current] = ""
        elif current is not None:
            records[current] += line.strip()
    return records


def test_parse_multidigit_allele_ids_safely():
    module = load_module()
    record = [
        "chrX", "100", ".", "A", "C,G,T,C,G,T,C,G,T,C", "99", "PASS", ".", "GT", "10/10"
    ]
    column = module.get_matrix_column(record, num_samples=1, resolve_IUPAC=False)
    assert column == "C"
    assert module.parse_genotype_alleles("10/10") == ["10", "10"]


def test_cli_audit_filters_and_outputs(tmp_path):
    outdir = tmp_path / "out"
    cmd = [
        sys.executable,
        str(SCRIPT),
        "-i", str(DATA / "synthetic_basic.vcf"),
        "--output-folder", str(outdir),
        "--output-prefix", "synthetic_basic",
        "--filter-pass-only",
        "--biallelic-only",
        "--exclude-bed", str(DATA / "exclude_regions.bed"),
        "--thin-distance", "25",
        "--audit-summary",
        "--write-used-sites",
        "-f",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    assert "SNPs that passed the filters: 2" in result.stdout

    prefix = outdir / "synthetic_basic.min4"
    fasta = read_fasta(prefix.with_suffix(".min4.fasta")) if False else read_fasta(str(prefix) + ".fasta")
    assert fasta["Sapria_A"] == "AT"
    assert fasta["Sapria_B"] == "RY"
    assert fasta["Rafflesia_A"] == "GC"
    assert fasta["Rafflesia_B"] == "AT"

    summary = json.loads(Path(str(prefix) + ".audit_summary.json").read_text())
    counts = summary["counts"]
    assert counts["records_processed"] == 6
    assert counts["records_accepted_nucleotide"] == 2
    assert counts["records_skipped_filter_not_pass"] == 1
    assert counts["records_skipped_not_biallelic"] == 1
    assert counts["records_skipped_exclude_bed"] == 1
    assert counts["records_skipped_thin_distance"] == 1

    used_sites = Path(str(prefix) + ".used_sites.tsv").read_text().splitlines()
    assert used_sites == ["#CHROM\tPOS\tNUM_SAMPLES", "chr1\t10\t4", "chr1\t40\t4"]


def test_seeded_iupac_resolution_is_reproducible(tmp_path):
    outputs = []
    for run_id in ("a", "b"):
        outdir = tmp_path / run_id
        cmd = [
            sys.executable,
            str(SCRIPT),
            "-i", str(DATA / "synthetic_basic.vcf"),
            "--output-folder", str(outdir),
            "--output-prefix", "seeded",
            "--filter-pass-only",
            "--biallelic-only",
            "--seed", "42",
            "--resolve-IUPAC",
            "-f",
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        outputs.append(Path(str(outdir / "seeded.min4") + ".fasta").read_text())
    assert outputs[0] == outputs[1]
