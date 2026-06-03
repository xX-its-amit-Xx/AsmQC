# AsmQC — GFF3 / BED export of flagged regions.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Turn the located QC features into a GFF3 and a BED track.

Features with genome coordinates are exported so a reviewer can inspect them in a
browser alongside the assembly:

* **contamination** — FCS-GX flagged spans (action, taxon).
* **telomere** — detected telomeric repeat arrays at scaffold ends.
* **gap** — assembly gaps (runs of N), recovered from the FASTA on demand.

GFF3 is 1-based inclusive; BED is 0-based half-open — both are emitted from the
same internal 1-based feature list.
"""
from __future__ import annotations

import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

from asmqc.fasta_stats import iter_fasta
from asmqc.models import QCReport

GFF_SOURCE = "AsmQC"


@dataclass
class TrackFeature:
    seq_id: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    ftype: str
    name: str
    score: float | None = None
    attributes: dict[str, str] = field(default_factory=dict)


def iter_gaps(fasta_path: str | Path, min_gap_len: int = 10) -> Iterator[tuple[str, int, int]]:
    """Yield ``(seq_id, start_1based, end_1based)`` for each assembly gap."""
    gap_re = re.compile(r"N{%d,}" % max(1, int(min_gap_len)))
    for name, seq in iter_fasta(fasta_path):
        upper = seq.upper()
        for m in gap_re.finditer(upper):
            yield name, m.start() + 1, m.end()


def build_features(
    report: QCReport,
    fasta_path: str | Path | None = None,
    min_gap_len: int = 10,
    include_gaps: bool = True,
) -> list[TrackFeature]:
    """Collect all coordinate-bearing QC features into one list."""
    features: list[TrackFeature] = []

    # Contamination (FCS provides coordinates; Kraken does not).
    if report.contamination and report.contamination.source == "fcs":
        for i, h in enumerate(report.contamination.hits):
            if h.start is None or h.end is None:
                continue
            features.append(TrackFeature(
                seq_id=h.seq_id,
                start=max(1, h.start),
                end=h.end,
                ftype="contamination",
                name=f"{h.action or 'contam'}_{i + 1}",
                score=h.coverage,
                attributes={k: v for k, v in {
                    "action": h.action or "",
                    "taxon": h.taxon or "",
                    "division": h.note or "",
                }.items() if v},
            ))

    # Telomeres.
    if report.telomere:
        win = report.telomere.window_size or 10000
        for s in report.telomere.scaffolds:
            length = s.length or 0
            if s.start_telomere:
                features.append(TrackFeature(
                    seq_id=s.name, start=1, end=min(win, length or win),
                    ftype="telomere", name=f"{s.name}_telomere_start",
                    score=float(s.start_repeats),
                    attributes={"end": "start", "repeats": str(s.start_repeats)},
                ))
            if s.end_telomere and length:
                features.append(TrackFeature(
                    seq_id=s.name, start=max(1, length - win + 1), end=length,
                    ftype="telomere", name=f"{s.name}_telomere_end",
                    score=float(s.end_repeats),
                    attributes={"end": "end", "repeats": str(s.end_repeats)},
                ))

    # Gaps (recovered from FASTA).
    if include_gaps and fasta_path and Path(fasta_path).exists():
        for j, (seq_id, start, end) in enumerate(iter_gaps(fasta_path, min_gap_len)):
            features.append(TrackFeature(
                seq_id=seq_id, start=start, end=end, ftype="gap",
                name=f"gap_{j + 1}", attributes={"length": str(end - start + 1)},
            ))

    features.sort(key=lambda f: (f.seq_id, f.start, f.end))
    return features


def _gff_escape(value: str) -> str:
    return (value.replace("%", "%25").replace(";", "%3B").replace("=", "%3D")
            .replace("&", "%26").replace(",", "%2C").replace("\t", "%09")
            .replace("\n", "%0A"))


# Characters allowed unescaped in a GFF3 seqid column (per the spec).
_SEQID_OK = re.compile(r"[A-Za-z0-9.:^*$@!+_?\-|]")


def _escape_seqid(value: str) -> str:
    """Percent-encode anything outside the GFF3-allowed seqid character set so a
    space/tab/newline in a sequence id cannot shift the tab-delimited columns."""
    return "".join(
        ch if _SEQID_OK.match(ch) else "%%%02X" % ord(ch)
        for ch in str(value)
    )


def _bed_token(value: str) -> str:
    """A BED field is tab-delimited and whitespace-free: collapse any whitespace
    (incl. embedded tabs/newlines) so a token cannot break the column layout."""
    return re.sub(r"\s+", "_", str(value).strip()) or "."


def write_gff3(features: list[TrackFeature], path: str | Path) -> Path:
    path = Path(path)
    lines = ["##gff-version 3"]
    for f in features:
        attrs = [f"ID={_gff_escape(f.name)}", f"Name={_gff_escape(f.name)}",
                 f"flag_type={_gff_escape(f.ftype)}"]
        for k, v in f.attributes.items():
            attrs.append(f"{_gff_escape(k)}={_gff_escape(str(v))}")
        score = "." if f.score is None else f"{f.score:g}"
        lines.append("\t".join([
            _escape_seqid(f.seq_id), GFF_SOURCE, f.ftype, str(f.start), str(f.end),
            score, ".", ".", ";".join(attrs),
        ]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_bed(features: list[TrackFeature], path: str | Path) -> Path:
    path = Path(path)
    lines = ['track name="AsmQC flags" description="AsmQC curation flags" itemRgb="On"']
    color = {"contamination": "214,39,40", "telomere": "51,160,44", "gap": "120,120,120"}
    for f in features:
        rgb = color.get(f.ftype, "31,120,180")
        score = 0 if f.score is None else max(0, min(1000, int(f.score)))
        # BED: 0-based start, half-open end.
        lines.append("\t".join([
            _bed_token(f.seq_id), str(f.start - 1), str(f.end),
            f"{f.ftype}:{_bed_token(f.name)}",
            str(score), ".", str(f.start - 1), str(f.end), rgb,
        ]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_tracks(
    report: QCReport,
    out_dir: str | Path,
    fasta_path: str | Path | None = None,
    min_gap_len: int = 10,
    include_gaps: bool = True,
) -> dict[str, Path]:
    """Write ``flags.gff3`` and ``flags.bed``; return the written paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    features = build_features(report, fasta_path, min_gap_len, include_gaps)
    return {
        "gff3": write_gff3(features, out_dir / "flags.gff3"),
        "bed": write_bed(features, out_dir / "flags.bed"),
        "n_features": len(features),  # type: ignore[dict-item]
    }
