# AsmQC — parsers for upstream QC tool outputs.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Parsers for each upstream tool's native output format.

Every ``parse_*`` function takes a path and returns either the populated result
model or ``None`` (when the file is missing/empty/unrecognisable).  They never
raise on malformed content — they collect what they can and leave the rest
``None`` — so the report degrades gracefully.
"""
from asmqc.parsers.busco import parse_busco
from asmqc.parsers.contamination import parse_contamination, parse_fcs, parse_kraken
from asmqc.parsers.gfastats import parse_gfastats
from asmqc.parsers.merqury import parse_merqury
from asmqc.parsers.quast import parse_quast
from asmqc.parsers.tidk import parse_tidk

__all__ = [
    "parse_busco",
    "parse_merqury",
    "parse_quast",
    "parse_gfastats",
    "parse_tidk",
    "parse_contamination",
    "parse_fcs",
    "parse_kraken",
]
