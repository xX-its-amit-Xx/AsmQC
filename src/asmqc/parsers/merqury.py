# AsmQC — Merqury output parser (.qv, .completeness.stats, .spectra-cn.hist).
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Parse Merqury k-mer evaluation outputs.

Formats verified against the Merqury source (github.com/marbl/merqury,
``eval/qv.sh`` and ``eval/spectra-cn.sh``) and the official wiki:

* ``<prefix>.qv`` — TAB, **no header**, 5 cols::

      col    682359    123511626    35.1183    0.000307729
      name   asm-only  total        QV         error-rate

* ``<prefix>.completeness.stats`` — TAB, **no header**, 5 cols::

      col    all    104975080    125303808    83.7764
      name   set    found        total        completeness%

* ``<prefix>.spectra-cn.hist`` — TAB, **with header** ``Copies\\tkmer_multiplicity\\tCount``.
  ``<prefix>.spectra-asm.hist`` is the same shape with header column ``Assembly``.

``parse_merqury`` accepts a directory (auto-discovers the files), a prefix, or any
one of the three files (it then looks for siblings sharing the prefix).
"""
from __future__ import annotations

from pathlib import Path

from asmqc.models import MerquryResult, SpectraPoint
from asmqc.parsers.common import open_text, to_float, to_int


def parse_qv_file(path: str | Path) -> tuple[float | None, dict[str, float], dict]:
    """Return ``(primary_qv, per_assembly_qv, extra)`` from a ``.qv`` file."""
    per_assembly: dict[str, float] = {}
    rows: list[tuple[str, int | None, int | None, float | None, float | None]] = []
    try:
        with open_text(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    parts = line.split()
                if len(parts) < 4:
                    continue
                name = parts[0].strip()
                if not name or name.lower() == "assembly":
                    continue
                asm_only = to_int(parts[1])
                total = to_int(parts[2])
                qv = to_float(parts[3])
                err = to_float(parts[4]) if len(parts) > 4 else None
                if qv is None:
                    continue
                rows.append((name, asm_only, total, qv, err))
                per_assembly[name] = qv
    except OSError:
        return None, {}, {}

    if not rows:
        return None, {}, {}
    # Prefer a combined 'Both' row, else the first assembly row.
    chosen = next((r for r in rows if r[0].lower() == "both"), rows[0])
    extra = {
        "kmers_asm_only": chosen[1],
        "kmers_total": chosen[2],
        "error_rate": chosen[4],
        "assembly_label": chosen[0],
    }
    return chosen[3], per_assembly, extra


def parse_completeness_file(path: str | Path) -> float | None:
    """Return the k-mer completeness percentage (region ``all`` preferred)."""
    best: float | None = None
    try:
        with open_text(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    parts = line.split()
                if len(parts) < 4:
                    continue
                name = parts[0].strip()
                if not name or name.lower() == "assembly":
                    continue
                region = parts[1].strip().lower()
                completeness = to_float(parts[4]) if len(parts) > 4 else None
                if completeness is None:
                    found = to_float(parts[2])
                    total = to_float(parts[3])
                    if found is not None and total:
                        completeness = 100.0 * found / total
                if completeness is None:
                    continue
                if region == "all":
                    return completeness
                if best is None:
                    best = completeness
    except OSError:
        return None
    return best


def parse_spectra_file(path: str | Path) -> list[SpectraPoint]:
    """Parse a ``.spectra-cn.hist`` / ``.spectra-asm.hist`` histogram."""
    points: list[SpectraPoint] = []
    try:
        with open_text(path) as fh:
            for i, line in enumerate(fh):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                copy_class = parts[0].strip()
                # Skip the header row (first column 'Copies' or 'Assembly').
                if i == 0 and copy_class.lower() in ("copies", "assembly"):
                    continue
                mult = to_int(parts[1])
                count = to_int(parts[2])
                if mult is None or count is None:
                    continue
                points.append(SpectraPoint(copy_class=copy_class, multiplicity=mult, count=count))
    except OSError:
        return []
    return points


def _discover(path: Path) -> dict[str, Path | None]:
    """Find the .qv / .completeness.stats / spectra histogram for a path."""
    found: dict[str, Path | None] = {"qv": None, "completeness": None, "spectra": None}

    def classify(p: Path) -> None:
        n = p.name.lower()
        if n.endswith(".qv") and found["qv"] is None:
            found["qv"] = p
        elif n.endswith(".completeness.stats") and found["completeness"] is None:
            found["completeness"] = p
        elif (n.endswith(".spectra-cn.hist") or n.endswith(".spectra-asm.hist")):
            # Prefer spectra-cn over spectra-asm.
            if found["spectra"] is None or n.endswith(".spectra-cn.hist"):
                found["spectra"] = p

    if path.is_dir():
        for p in sorted(path.iterdir()):
            if p.is_file():
                classify(p)
    elif path.is_file():
        classify(path)
        # Look for siblings sharing the prefix (strip the known suffix).
        name = path.name
        prefix = name
        for suf in (".qv", ".completeness.stats", ".spectra-cn.hist", ".spectra-asm.hist"):
            if name.endswith(suf):
                prefix = name[: -len(suf)]
                break
        for p in sorted(path.parent.iterdir()):
            if p.is_file() and p.name.startswith(prefix):
                classify(p)
    else:
        # Treat as a prefix path that does not itself exist.
        parent = path.parent
        prefix = path.name
        if parent.is_dir():
            for p in sorted(parent.iterdir()):
                if p.is_file() and p.name.startswith(prefix):
                    classify(p)
    return found


def parse_merqury(
    path: str | Path | None = None,
    *,
    qv: str | Path | None = None,
    completeness: str | Path | None = None,
    spectra: str | Path | None = None,
    kmer_size: int | None = None,
) -> MerquryResult | None:
    """Parse Merqury outputs from a directory, prefix, or explicit files."""
    files: dict[str, Path | None] = {"qv": None, "completeness": None, "spectra": None}
    if path is not None:
        files = _discover(Path(path))
    if qv:
        files["qv"] = Path(qv)
    if completeness:
        files["completeness"] = Path(completeness)
    if spectra:
        files["spectra"] = Path(spectra)

    if not any(files.values()):
        return None

    result = MerquryResult(kmer_size=kmer_size)
    if files["qv"]:
        primary, per_assembly, extra = parse_qv_file(files["qv"])
        result.qv = primary
        result.per_assembly_qv = per_assembly
        result.kmers_asm_only = extra.get("kmers_asm_only")
        result.kmers_total = extra.get("kmers_total")
        result.error_rate = extra.get("error_rate")
        result.assembly_label = extra.get("assembly_label")
    if files["completeness"]:
        result.completeness_pct = parse_completeness_file(files["completeness"])
    if files["spectra"]:
        result.spectra = parse_spectra_file(files["spectra"])

    has_data = (
        result.qv is not None
        or result.completeness_pct is not None
        or bool(result.spectra)
    )
    return result if has_data else None
