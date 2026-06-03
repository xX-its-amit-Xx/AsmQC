# AsmQC — matplotlib figures embedded as base64 PNGs.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Generate the report figures and return them as base64-encoded PNGs.

Every builder is defensive: if the required data is absent or plotting fails it
returns ``None`` so the report still renders.  Figures are produced with the
non-interactive Agg backend, so no display is required (CI-safe).
"""
from __future__ import annotations

import base64
import io

import matplotlib

matplotlib.use("Agg")  # headless backend — must precede pyplot import
import matplotlib.pyplot as plt  # noqa: E402

from asmqc.config import AsmQCConfig  # noqa: E402
from asmqc.models import (  # noqa: E402
    AssemblyStats,
    BuscoResult,
    MerquryResult,
    QCReport,
    TelomereResult,
)

# Consistent palette.
_BUSCO_COLORS = {"S": "#1f78b4", "D": "#a6cee3", "F": "#fdbf6f", "M": "#e31a1c"}
_SPECTRA_COLORS = {
    "read-only": "#000000", "0": "#000000", "absent": "#000000",
    "1": "#e31a1c", "2": "#1f78b4", "3": "#33a02c", "4": "#6a3d9a", ">4": "#ff7f00",
}


def _fig_to_base64(fig, dpi: int = 110) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def build_plots(report: QCReport, config: AsmQCConfig) -> dict[str, str]:
    """Return ``{plot_name: base64_png}`` for all plots that could be built."""
    dpi = int(config.get("report", "figure_dpi", 110))
    max_seqs = int(config.get("report", "max_seqs_in_plots", 200))
    plots: dict[str, str | None] = {}
    plots["spectra_cn"] = spectra_cn_plot(report.merqury, dpi) if report.merqury else None
    plots["cumulative_length"] = cumulative_length_plot(report.stats, dpi) if report.stats else None
    plots["busco"] = busco_plot(report.busco, dpi) if report.busco else None
    chrom_min = int(config.get("telomere", "chromosome_min_length_bp", 1_000_000))
    plots["telomere_map"] = (
        telomere_map_plot(report.telomere, dpi, max_seqs, chrom_min)
        if report.telomere else None
    )
    return {k: v for k, v in plots.items() if v}


def spectra_cn_plot(m: MerquryResult, dpi: int = 110) -> str | None:
    if not m or not m.spectra:
        return None
    try:
        by_class: dict[str, list[tuple[int, int]]] = {}
        for pt in m.spectra:
            by_class.setdefault(pt.copy_class, []).append((pt.multiplicity, pt.count))
        if not by_class:
            return None
        fig, ax = plt.subplots(figsize=(7, 4.2))
        order = ["read-only", "0", "absent", "1", "2", "3", "4", ">4"]
        classes = sorted(by_class, key=lambda c: order.index(c) if c in order else 99)
        max_mult = 1
        for cls in classes:
            pts = sorted(by_class[cls])
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            if xs:
                max_mult = max(max_mult, max(xs))
            color = _SPECTRA_COLORS.get(cls, None)
            ax.fill_between(xs, ys, alpha=0.35, color=color)
            ax.plot(xs, ys, label=f"copies: {cls}", color=color, linewidth=1.2)
        # Trim the long high-multiplicity tail for readability.
        ax.set_xlim(0, min(max_mult, _percentile_cap(by_class)))
        ax.set_xlabel("k-mer multiplicity in reads")
        ax.set_ylabel("number of k-mers")
        ax.set_title("Merqury copy-number (spectra-cn) spectrum")
        ax.legend(fontsize=8, ncol=2)
        ax.margins(y=0.02)
        return _fig_to_base64(fig, dpi)
    except Exception:
        return None


def _percentile_cap(by_class: dict[str, list[tuple[int, int]]]) -> int:
    # Cap x-axis a bit past where the bulk of k-mers live.
    mults = [mult for pts in by_class.values() for mult, count in pts if count > 0]
    if not mults:
        return 1
    mults.sort()
    return max(mults[min(len(mults) - 1, int(len(mults) * 0.98))], 10)


def cumulative_length_plot(s: AssemblyStats, dpi: int = 110) -> str | None:
    lengths = s.sorted_lengths if s and s.sorted_lengths else None
    if not lengths:
        return None
    try:
        total = sum(lengths)
        if total <= 0:
            return None
        # Nx curve: x = cumulative % of genome, y = scaffold length (log).
        xs = [0.0]
        ys = [lengths[0]]
        cum = 0
        for length in lengths:
            cum += length
            xs.append(100.0 * cum / total)
            ys.append(length)
        fig, ax = plt.subplots(figsize=(7, 4.2))
        ax.step(xs, ys, where="post", color="#1f78b4", linewidth=1.6)
        ax.set_yscale("log")
        ax.set_xlabel("cumulative % of assembly length")
        ax.set_ylabel("scaffold length (bp, log scale)")
        ax.set_title("Cumulative scaffold-length (Nx) curve")
        if s.n50:
            ax.axvline(50, color="grey", linestyle="--", linewidth=1)
            ax.axhline(s.n50, color="#e31a1c", linestyle=":", linewidth=1)
            ax.annotate(
                f"N50 = {_human_bp(s.n50)}", xy=(50, s.n50),
                xytext=(52, s.n50 * 1.3), fontsize=8, color="#e31a1c",
            )
        ax.grid(True, which="both", alpha=0.25)
        return _fig_to_base64(fig, dpi)
    except Exception:
        return None


def busco_plot(b: BuscoResult, dpi: int = 110) -> str | None:
    if not b:
        return None
    s = b.single_pct
    d = b.duplicated_pct
    f = b.fragmented_pct
    m = b.missing_pct
    if all(v is None for v in (s, d, f, m)):
        return None
    s, d, f, m = (v or 0.0 for v in (s, d, f, m))
    try:
        fig, ax = plt.subplots(figsize=(7, 1.9))
        left = 0.0
        for value, key, label in (
            (s, "S", "Complete single-copy (S)"),
            (d, "D", "Complete duplicated (D)"),
            (f, "F", "Fragmented (F)"),
            (m, "M", "Missing (M)"),
        ):
            ax.barh(0, value, left=left, color=_BUSCO_COLORS[key],
                    edgecolor="white", label=f"{label}: {value:.1f}%")
            if value >= 4:
                ax.text(left + value / 2, 0, f"{value:.1f}", va="center", ha="center",
                        fontsize=8, color="white")
            left += value
        ax.set_xlim(0, max(100.0, left))
        ax.set_yticks([])
        ax.set_xlabel("% of BUSCO groups")
        title = "BUSCO completeness"
        if b.lineage:
            title += f" ({b.lineage})"
        ax.set_title(title)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.45), ncol=2, fontsize=8,
                  frameon=False)
        return _fig_to_base64(fig, dpi)
    except Exception:
        return None


def telomere_map_plot(
    t: TelomereResult, dpi: int = 110, max_seqs: int = 200, chrom_min: int = 1_000_000
) -> str | None:
    if not t or not t.scaffolds:
        return None
    try:
        chroms = [s for s in t.scaffolds if (s.length or 0) >= chrom_min] or t.scaffolds
        chroms = sorted(chroms, key=lambda s: -(s.length or 0))[: min(max_seqs, 50)]
        if not chroms:
            return None
        n = len(chroms)
        max_len = max((s.length or 1) for s in chroms)
        fig, ax = plt.subplots(figsize=(7, max(2.2, 0.28 * n + 1)))
        for i, s in enumerate(reversed(chroms)):
            length = s.length or max_len
            width = length / max_len
            ax.barh(i, width, height=0.55, color="#dddddd", edgecolor="#999999", zorder=1)
            cap = 0.018
            start_color = "#33a02c" if s.start_telomere else "#e31a1c"
            end_color = "#33a02c" if s.end_telomere else "#e31a1c"
            ax.barh(i, cap, left=0, height=0.55, color=start_color, zorder=2)
            ax.barh(i, cap, left=width - cap, height=0.55, color=end_color, zorder=2)
        ax.set_yticks(range(n))
        ax.set_yticklabels([s.name for s in reversed(chroms)], fontsize=7)
        ax.set_xlim(0, 1.02)
        ax.set_xlabel("scaffold (length-normalised); green = telomere found, red = none")
        title = "Per-scaffold telomere map"
        if t.repeat_motif:
            title += f" (motif {t.repeat_motif})"
        ax.set_title(title)
        ax.spines[["top", "right"]].set_visible(False)
        return _fig_to_base64(fig, dpi)
    except Exception:
        return None


def _human_bp(value: float) -> str:
    if value >= 1e9:
        return f"{value / 1e9:.2f} Gb"
    if value >= 1e6:
        return f"{value / 1e6:.2f} Mb"
    if value >= 1e3:
        return f"{value / 1e3:.1f} kb"
    return f"{int(value)} bp"
