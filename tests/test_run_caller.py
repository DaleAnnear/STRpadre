from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import run_caller


def longtr_command(options: dict) -> list[str]:
    config = {
        "threads": 1,
        "phased_reads_expected": False,
        "additional_args": [],
        "options": options,
    }
    args = Namespace(
        alignment=Path("reads.bam"),
        reference=Path("reference.fa"),
        catalog=Path("catalog.bed"),
        sample_id="SAMPLE",
    )
    command, _, _ = run_caller.build_command("longtr", config, args, Path("output"))
    return command


def test_longtr_uses_sample_as_library_by_default() -> None:
    assert "--lib-from-samp" in longtr_command({})


def test_longtr_can_use_read_group_library_tags() -> None:
    assert "--lib-from-samp" not in longtr_command({"use_lb_tags": True})
