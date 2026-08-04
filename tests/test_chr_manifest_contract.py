from __future__ import annotations

import importlib.util
from pathlib import Path


def load_adapter():
    path = Path(__file__).parents[1] / "bin" / "catalog_adapter.py"
    spec = importlib.util.spec_from_file_location("catalog_adapter", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_catalog_adapter_requires_chr_not_legacy_contig(tmp_path: Path, capsys) -> None:
    adapter = load_adapter()
    manifest = tmp_path / "legacy.tsv"
    fai = tmp_path / "reference.fa.fai"
    manifest.write_text(
        "locus_id\tcontig\tstart\tend\tcoordinate_system\treference_build\tprimary_motif\n"
        "L1\tchr1\t10\t16\t0-based-half-open\ttest-build\tCAG\n"
    )
    fai.write_text("chr1\t100\t0\t0\t0\n")

    assert adapter.main([
        "--locus-manifest", str(manifest),
        "--reference-fai", str(fai),
        "--callers", "longtr",
        "--output-dir", str(tmp_path / "catalogs"),
    ]) == 2
    assert "chr" in capsys.readouterr().err
