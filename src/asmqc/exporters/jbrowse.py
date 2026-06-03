# AsmQC — JBrowse2 config exporter.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Write a minimal JBrowse2 ``config.json`` so a reviewer can load the assembly
with the AsmQC flag track overlaid.

The config references the assembly FASTA via an ``IndexedFastaAdapter`` (which
needs a ``.fai`` index — run ``samtools faidx assembly.fa``) and the AsmQC
``flags.gff3`` via a plain ``Gff3Adapter``.  URIs are relative filenames so the
config can sit in the same directory as the assembly and track.
"""
from __future__ import annotations

import json
import re
from pathlib import Path


def _sanitise(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "assembly"


def build_jbrowse_config(
    assembly_name: str,
    fasta_uri: str,
    gff3_uri: str = "flags.gff3",
    fai_uri: str | None = None,
    bed_uri: str | None = None,
) -> dict:
    asm_id = _sanitise(assembly_name)
    fai_uri = fai_uri or f"{fasta_uri}.fai"

    config: dict = {
        "assemblies": [
            {
                "name": asm_id,
                "sequence": {
                    "type": "ReferenceSequenceTrack",
                    "trackId": f"{asm_id}-ref",
                    "adapter": {
                        "type": "IndexedFastaAdapter",
                        "fastaLocation": {"uri": fasta_uri, "locationType": "UriLocation"},
                        "faiLocation": {"uri": fai_uri, "locationType": "UriLocation"},
                    },
                },
            }
        ],
        "tracks": [
            {
                "type": "FeatureTrack",
                "trackId": "asmqc-flags",
                "name": "AsmQC curation flags",
                "category": ["AsmQC"],
                "assemblyNames": [asm_id],
                "adapter": {
                    "type": "Gff3Adapter",
                    "gffLocation": {"uri": gff3_uri, "locationType": "UriLocation"},
                },
            }
        ],
        "defaultSession": {
            "name": "AsmQC review",
            "views": [
                {
                    "type": "LinearGenomeView",
                    "tracks": [
                        {
                            "type": "FeatureTrack",
                            "configuration": "asmqc-flags",
                            "displays": [
                                {"type": "LinearBasicDisplay", "configuration": "asmqc-flags-LinearBasicDisplay"}
                            ],
                        }
                    ],
                }
            ],
        },
    }
    if bed_uri:
        config["tracks"].append({
            "type": "FeatureTrack",
            "trackId": "asmqc-flags-bed",
            "name": "AsmQC curation flags (BED)",
            "category": ["AsmQC"],
            "assemblyNames": [asm_id],
            "adapter": {
                "type": "BedAdapter",
                "bedLocation": {"uri": bed_uri, "locationType": "UriLocation"},
            },
        })
    return config


def write_jbrowse_config(
    out_dir: str | Path,
    assembly_name: str,
    fasta_path: str | Path,
    gff3_name: str = "flags.gff3",
    bed_name: str | None = None,
) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Reference the FASTA by basename so the config is portable next to it.
    fasta_uri = Path(fasta_path).name
    config = build_jbrowse_config(
        assembly_name=assembly_name,
        fasta_uri=fasta_uri,
        gff3_uri=gff3_name,
        bed_uri=bed_name,
    )
    path = out_dir / "jbrowse2_config.json"
    path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    return path
