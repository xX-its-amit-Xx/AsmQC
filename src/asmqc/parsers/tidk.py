# AsmQC — tidk telomere window parser.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Parse tidk ``*_telomeric_repeat_windows.tsv`` (and the bedgraph variant).

Format verified against tidk source (``tolkit/telomeric-identifier``,
``src/search.rs`` / ``src/finder.rs``):

* TSV, **with header** ``id\\twindow\\tforward_repeat_number\\treverse_repeat_number\\ttelomeric_repeat``.
* ``window`` is the **end coordinate** of the window (``window_size``, 2*size, ...,
  last clamped to the scaffold length) — so ``max(window)`` per id ≈ scaffold length.
* bedgraph variant: **no header**, 4 cols ``id\\tstart\\tend\\t(fwd+rev)``.

Telomere inference (biological convention): a telomere sits at a chromosome end,
so we look only at the first/last ``terminal_windows`` of each scaffold.  A
terminal window whose motif count reaches ``min_window_count`` marks a telomere
there.  (At a scaffold *start* the leading-strand telomere shows up as the
reverse-complement motif; at the *end* as the forward motif — we take the max of
the two to stay orientation-robust.)
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from asmqc.models import TelomereResult, TelomereScaffold
from asmqc.parsers.common import open_text, to_int

_TSV_HEADER = ("id", "window", "forward_repeat_number", "reverse_repeat_number", "telomeric_repeat")


def parse_tidk(
    path: str | Path,
    *,
    min_window_count: int = 50,
    terminal_windows: int = 2,
    chromosome_min_length: int = 0,
) -> TelomereResult | None:
    p = Path(path)
    if not p.exists():
        return None

    # rows[id] = list of (window_end, fwd, rev)
    rows: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    motif: str | None = None
    is_bedgraph = p.name.lower().endswith(".bedgraph")

    try:
        with open_text(p) as fh:
            for i, line in enumerate(fh):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                if i == 0 and parts[0].strip().lower() == "id" and not is_bedgraph:
                    continue  # header
                seq_id = parts[0].strip()
                if not seq_id:
                    continue
                if is_bedgraph:
                    # id, start, end, summed count
                    window_end = to_int(parts[2])
                    total = to_int(parts[3]) or 0
                    fwd, rev = total, 0
                else:
                    window_end = to_int(parts[1])
                    fwd = to_int(parts[2]) or 0
                    rev = to_int(parts[3]) or 0
                if window_end is None:
                    continue  # header or junk line — skip before reading motif
                # Capture the motif only from a validated data row, so a header
                # not on line 0 (e.g. after a leading blank line) can't poison it.
                if not is_bedgraph and len(parts) >= 5 and motif is None:
                    motif = parts[4].strip() or None
                rows[seq_id].append((window_end, fwd, rev))
    except OSError:
        return None

    if not rows:
        return None

    window_size = _infer_window_size(rows)
    scaffolds: list[TelomereScaffold] = []
    for seq_id, windows in rows.items():
        windows.sort(key=lambda w: w[0])
        length = windows[-1][0] if windows else None
        k = max(1, terminal_windows)
        nwin = len(windows)
        if nwin == 1:
            # Can't separate the two ends from one window: use the strand
            # convention (reverse-complement motif marks the 5' start, forward
            # motif the 3' end) so a single-end telomere isn't counted twice.
            _, fwd0, rev0 = windows[0]
            start_count, end_count = rev0, fwd0
        else:
            # Keep head and tail disjoint so a one-end telomere on a scaffold
            # with few windows isn't double-counted. Take max(fwd,rev) per
            # terminal window for orientation robustness.
            kk = min(k, nwin // 2)
            head = windows[:kk]
            tail = windows[-kk:]
            start_count = max((max(f, r) for _, f, r in head), default=0)
            end_count = max((max(f, r) for _, f, r in tail), default=0)
        total_repeats = sum(f + r for _, f, r in windows)
        scaffolds.append(
            TelomereScaffold(
                name=seq_id,
                length=length,
                start_telomere=start_count >= min_window_count,
                end_telomere=end_count >= min_window_count,
                start_repeats=start_count,
                end_repeats=end_count,
                total_repeats=total_repeats,
            )
        )

    scaffolds.sort(key=lambda s: (-(s.length or 0), s.name))
    return TelomereResult(
        repeat_motif=motif,
        window_size=window_size,
        scaffolds=scaffolds,
        chromosome_min_length=chromosome_min_length,
    )


def _infer_window_size(rows: dict[str, list[tuple[int, int, int]]]) -> int | None:
    from collections import Counter

    diffs: Counter = Counter()
    for windows in rows.values():
        ws = sorted(w[0] for w in windows)
        for a, b in zip(ws, ws[1:]):
            if b - a > 0:
                diffs[b - a] += 1
    if diffs:
        return diffs.most_common(1)[0][0]
    # Single-window scaffolds: the window value itself is the size (or length).
    for windows in rows.values():
        if windows:
            return windows[0][0]
    return None
