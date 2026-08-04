from __future__ import annotations

import csv
import gzip
import importlib.util
from pathlib import Path


def load_adapter():
    path = Path(__file__).parents[1] / "bin" / "catalog_adapter.py"
    spec = importlib.util.spec_from_file_location("catalog_adapter", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FakePysam:
    @staticmethod
    def tabix_compress(source: str, destination: str, force: bool = False) -> None:
        del force
        with open(source, "rb") as input_handle, gzip.open(destination, "wb") as output_handle:
            output_handle.write(input_handle.read())

    @staticmethod
    def tabix_index(path: str, preset: str, force: bool = False) -> None:
        assert preset == "bed"
        assert force
        Path(f"{path}.tbi").write_bytes(b"test-tabix-index\n")


def write_manifest(path: Path, repeat_structure: str = "") -> None:
    path.write_text(
        "locus_id\tchr\tstart\tend\tcoordinate_system\treference_build\tprimary_motif\trepeat_structure\ttrgt_id\tlongtr_id\tatarva_id\tstrdust_id\n"
        f"LOC_B\tchr2\t4\t10\t0-based-half-open\ttest-build\tGAA\t{repeat_structure}\tTRGT_B\tLONG_B\tAT_B\tSTR_B\n"
        "LOC_A\tchr1\t10\t16\t0-based-half-open\ttest-build\tCAG\t\tTRGT_A\tLONG_A\tAT_A\tSTR_A\n"
    )


def test_adapter_emits_all_native_catalogs_and_provenance(tmp_path: Path, monkeypatch) -> None:
    adapter = load_adapter()
    monkeypatch.setattr(adapter, "pysam", FakePysam)
    manifest = tmp_path / "loci.tsv"
    fai = tmp_path / "reference.fa.fai"
    write_manifest(manifest)
    fai.write_text("chr1\t100\t0\t0\t0\nchr2\t100\t0\t0\t0\n")
    output = tmp_path / "catalogs"

    assert (
        adapter.main(
            [
                "--locus-manifest",
                str(manifest),
                "--reference-fai",
                str(fai),
                "--callers",
                "all",
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )

    assert (output / "trgt.bed").read_text().splitlines() == [
        "chr1\t10\t16\tID=TRGT_A;MOTIFS=CAG;STRUC=<TR>",
        "chr2\t4\t10\tID=TRGT_B;MOTIFS=GAA;STRUC=<TR>",
    ]
    assert (output / "longtr.bed").read_text().splitlines() == [
        "chr1\t11\t16\tCAG\tLONG_A",
        "chr2\t5\t10\tGAA\tLONG_B",
    ]
    assert (output / "strdust.bed").read_text().splitlines() == ["chr1\t10\t16", "chr2\t4\t10"]
    with gzip.open(output / "atarva.bed.gz", "rt") as handle:
        assert handle.read().splitlines() == [
            "#CHROM\tSTART\tEND\tMOTIF\tMOTIF_LEN",
            "chr1\t10\t16\tCAG\t3",
            "chr2\t4\t10\tGAA\t3",
        ]
    assert (output / "atarva.bed.gz.tbi").is_file()
    with (output / "catalog_locus_mapping.tsv").open(newline="") as handle:
        mappings = list(csv.DictReader(handle, delimiter="\t"))
    assert len(mappings) == 8
    assert {row["caller"] for row in mappings} == {"trgt", "longtr", "atarva", "strdust"}
    assert {row["catalog_coordinate_system"] for row in mappings if row["caller"] == "longtr"} == {
        "1-based-closed"
    }
    assert (
        '"reference_build": "test-build"' in (output / "catalog_adaptation_report.json").read_text()
    )


def test_adapter_refuses_complex_trgt_without_partial_output(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    adapter = load_adapter()
    monkeypatch.setattr(adapter, "pysam", FakePysam)
    manifest = tmp_path / "complex.tsv"
    fai = tmp_path / "reference.fa.fai"
    write_manifest(manifest, repeat_structure="(GAA)*CCG")
    fai.write_text("chr1\t100\t0\t0\t0\nchr2\t100\t0\t0\t0\n")
    output = tmp_path / "catalogs"

    assert (
        adapter.main(
            [
                "--locus-manifest",
                str(manifest),
                "--reference-fai",
                str(fai),
                "--callers",
                "trgt,longtr",
                "--output-dir",
                str(output),
            ]
        )
        == 2
    )
    assert not output.exists()
    assert "complex locus LOC_B" in capsys.readouterr().err


def test_adapter_rejects_duplicate_callers_before_writing(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    adapter = load_adapter()
    monkeypatch.setattr(adapter, "pysam", FakePysam)
    manifest = tmp_path / "loci.tsv"
    fai = tmp_path / "reference.fa.fai"
    write_manifest(manifest)
    fai.write_text("chr1\t100\t0\t0\t0\nchr2\t100\t0\t0\t0\n")
    output = tmp_path / "catalogs"

    assert (
        adapter.main(
            [
                "--locus-manifest",
                str(manifest),
                "--reference-fai",
                str(fai),
                "--callers",
                "trgt,trgt",
                "--output-dir",
                str(output),
            ]
        )
        == 2
    )
    assert not output.exists()
    assert "Duplicate caller names" in capsys.readouterr().err
