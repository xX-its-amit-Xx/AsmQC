# AsmQC — report writers (HTML, Markdown, JSON).
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Render a :class:`~asmqc.models.QCReport` to a self-contained HTML page, a
Markdown document, and a machine-readable JSON summary.

The HTML embeds CSS and base64 PNG figures inline, so ``report.html`` is a single
portable file with no external assets.
"""
from __future__ import annotations

import html
import json
from pathlib import Path

from asmqc.models import QCReport, Severity

_SEVERITY_COLOR = {
    "PASS": "#2e8b57",
    "INFO": "#6c757d",
    "NOTE": "#1f78b4",
    "WARN": "#e08e0b",
    "FLAG": "#d62728",
    "FAIL": "#7a1416",
}
_PLOT_TITLES = {
    "busco": "BUSCO completeness",
    "spectra_cn": "Merqury k-mer copy-number spectrum",
    "cumulative_length": "Cumulative scaffold-length (Nx) curve",
    "telomere_map": "Per-scaffold telomere map",
}


# --------------------------------------------------------------------------- #
# Value formatting
# --------------------------------------------------------------------------- #
def _bp(value: float | None) -> str:
    if value is None:
        return "—"
    value = float(value)
    if value >= 1e9:
        return f"{value / 1e9:.2f} Gb"
    if value >= 1e6:
        return f"{value / 1e6:.2f} Mb"
    if value >= 1e3:
        return f"{value / 1e3:.1f} kb"
    return f"{int(value):,} bp"


def _num(value: float | None, suffix: str = "", pct: bool = False) -> str:
    if value is None:
        return "—"
    if pct:
        return f"{value:.1f}%"
    if isinstance(value, float) and not value.is_integer():
        return f"{value:,.2f}{suffix}"
    return f"{int(value):,}{suffix}"


def _summary_rows(report: QCReport) -> list[tuple[str, str, str]]:
    """Return (group, metric, value) rows for the summary table."""
    rows: list[tuple[str, str, str]] = []
    s = report.stats
    if s:
        rows += [
            ("Contiguity", "Total length", _bp(s.total_length)),
            ("Contiguity", "Sequences (scaffolds)", _num(s.num_sequences)),
            ("Contiguity", "Contigs", _num(s.num_contigs)),
            ("Contiguity", "Scaffold N50", _bp(s.n50)),
            ("Contiguity", "Scaffold N90", _bp(s.n90)),
            ("Contiguity", "L50", _num(s.l50)),
            ("Contiguity", "Contig N50", _bp(s.contig_n50)),
            ("Contiguity", "Largest", _bp(s.largest)),
            ("Contiguity", "GC content", _num(s.gc_percent, pct=True)),
            ("Contiguity", "Gaps", _num(s.gap_count)),
            ("Contiguity", "Gap bases", _bp(s.gap_bases)),
            ("Contiguity", "N per 100 kbp", _num(s.n_per_100kbp)),
        ]
    b = report.busco
    if b:
        rows += [
            ("Completeness (BUSCO)", "Complete", _num(b.complete_pct, pct=True)),
            ("Completeness (BUSCO)", "Single-copy", _num(b.single_pct, pct=True)),
            ("Completeness (BUSCO)", "Duplicated", _num(b.duplicated_pct, pct=True)),
            ("Completeness (BUSCO)", "Fragmented", _num(b.fragmented_pct, pct=True)),
            ("Completeness (BUSCO)", "Missing", _num(b.missing_pct, pct=True)),
            ("Completeness (BUSCO)", "Lineage", b.lineage or "—"),
        ]
    m = report.merqury
    if m:
        rows += [
            ("Correctness (Merqury)", "Consensus QV", _num(m.qv)),
            ("Correctness (Merqury)", "k-mer completeness", _num(m.completeness_pct, pct=True)),
        ]
    t = report.telomere
    if t:
        rows += [
            ("Telomeres", "Scaffolds w/ both ends", _num(t.n_with_both)),
            ("Telomeres", "Scaffolds w/ any end", _num(t.n_with_any)),
            ("Telomeres", "Repeat motif", t.repeat_motif or "—"),
        ]
    c = report.contamination
    if c:
        rows += [
            ("Contamination", "Source", c.source.upper()),
            ("Contamination", "Sequences flagged", _num(c.n_sequences_flagged)),
            ("Contamination", "Flagged bases", _bp(c.total_flagged_bases)),
        ]
    return rows


# --------------------------------------------------------------------------- #
# JSON
# --------------------------------------------------------------------------- #
def write_json(report: QCReport, path: str | Path, include_per_seq: bool = False) -> Path:
    path = Path(path)
    path.write_text(
        json.dumps(report.to_dict(include_per_seq=include_per_seq), indent=2),
        encoding="utf-8",
    )
    return path


# --------------------------------------------------------------------------- #
# Markdown
# --------------------------------------------------------------------------- #
def _md(value: object) -> str:
    """Escape Markdown table/HTML-significant characters in an untrusted value.

    Prevents a stray ``|`` from corrupting a table row and neutralises inline
    HTML so the Markdown is safe to render to HTML downstream. (Flag prose is
    already sanitised at the source in ``flags.py``.)
    """
    return (str(value).replace("\\", "\\\\").replace("|", "\\|")
            .replace("<", "&lt;").replace(">", "&gt;")
            .replace("\r", " ").replace("\n", " "))


def render_markdown(report: QCReport) -> str:
    out: list[str] = []
    status = report.overall_severity.name
    out.append(f"# AsmQC report — {_md(report.assembly_name)}\n")
    out.append(f"**Overall status:** {status}  ")
    if report.generated_at:
        out.append(f"\n_Generated {report.generated_at} by AsmQC {report.asmqc_version}_\n")

    counts = report.counts_by_severity()
    badge = " · ".join(f"{k}: {v}" for k, v in counts.items() if v)
    if badge:
        out.append(f"\n**Flag summary:** {badge}\n")

    # Flags
    out.append("\n## Curation flags\n")
    if not report.flags:
        out.append("_No flags — no QC inputs evaluated._\n")
    for f in report.flags:
        out.append(f"### [{f.severity.name}] {f.title}\n")
        out.append(f"- **Observation:** {f.message}")
        if f.explanation:
            out.append(f"- **Likely cause:** {f.explanation}")
        if f.action:
            out.append(f"- **Suggested action:** {f.action}")
        out.append("")

    # Summary metrics
    out.append("\n## Metrics\n")
    out.append("| Group | Metric | Value |")
    out.append("| --- | --- | --- |")
    for group, metric, value in _summary_rows(report):
        out.append(f"| {_md(group)} | {_md(metric)} | {_md(value)} |")

    # Inputs
    if report.inputs:
        out.append("\n## Inputs\n")
        for key, value in report.inputs.items():
            out.append(f"- **{_md(key)}:** {_md(value)}")

    if report.warnings:
        out.append("\n## Parse warnings\n")
        for w in report.warnings:
            out.append(f"- {_md(w)}")

    out.append("\n---\n")
    out.append("Thresholds follow the Earth BioGenome Project / VGP 3C framework "
               "(Contiguity, Completeness, Correctness). See the AsmQC README for "
               "what each metric means and what good vs bad looks like.\n")
    return "\n".join(out) + "\n"


def write_markdown(report: QCReport, path: str | Path) -> Path:
    path = Path(path)
    path.write_text(render_markdown(report), encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #
def render_html(report: QCReport, plots: dict[str, str] | None = None) -> str:
    plots = plots or {}
    e = html.escape
    status = report.overall_severity.name
    status_color = _SEVERITY_COLOR.get(status, "#6c757d")

    parts: list[str] = []
    parts.append("<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>")
    parts.append("<meta name='viewport' content='width=device-width, initial-scale=1'>")
    parts.append(f"<title>AsmQC report — {e(report.assembly_name)}</title>")
    parts.append(f"<style>{_CSS}</style></head><body><div class='wrap'>")

    # Header
    parts.append("<header>")
    parts.append("<h1>AsmQC report</h1>")
    parts.append(f"<div class='subtitle'>{e(report.assembly_name)}</div>")
    parts.append(
        f"<div class='status' style='background:{status_color}'>Overall: {status}</div>"
    )
    meta = []
    if report.generated_at:
        meta.append(f"Generated {e(report.generated_at)}")
    if report.asmqc_version:
        meta.append(f"AsmQC {e(report.asmqc_version)}")
    if meta:
        parts.append(f"<div class='meta'>{' · '.join(meta)}</div>")
    # Flag count chips
    chips = []
    for sev in Severity:
        n = report.counts_by_severity().get(sev.name, 0)
        if n:
            color = _SEVERITY_COLOR[sev.name]
            chips.append(f"<span class='chip' style='background:{color}'>{sev.name}: {n}</span>")
    if chips:
        parts.append(f"<div class='chips'>{''.join(chips)}</div>")
    parts.append("</header>")

    # Flags
    parts.append("<section><h2>Curation flags</h2>")
    if not report.flags:
        parts.append("<p class='muted'>No flags — no QC inputs were evaluated.</p>")
    for f in report.flags:
        color = _SEVERITY_COLOR.get(f.severity.name, "#6c757d")
        parts.append("<div class='flag'>")
        parts.append(
            f"<div class='flag-head'><span class='sev' style='background:{color}'>"
            f"{f.severity.name}</span><span class='flag-title'>{e(f.title)}</span>"
            f"<span class='flag-cat'>{e(f.category)}</span></div>"
        )
        parts.append(f"<div class='flag-msg'>{e(f.message)}</div>")
        if f.explanation:
            parts.append(f"<div class='flag-why'><b>Likely cause.</b> {e(f.explanation)}</div>")
        if f.action:
            parts.append(f"<div class='flag-act'><b>Action.</b> {e(f.action)}</div>")
        parts.append("</div>")
    parts.append("</section>")

    # Plots
    if plots:
        parts.append("<section><h2>Figures</h2><div class='plots'>")
        for key in ("busco", "spectra_cn", "cumulative_length", "telomere_map"):
            if key in plots:
                title = _PLOT_TITLES.get(key, key)
                parts.append(
                    f"<figure><img alt='{e(title)}' "
                    f"src='data:image/png;base64,{plots[key]}'>"
                    f"<figcaption>{e(title)}</figcaption></figure>"
                )
        parts.append("</div></section>")

    # Metrics table
    parts.append("<section><h2>Metrics</h2><table class='metrics'>")
    parts.append("<thead><tr><th>Group</th><th>Metric</th><th>Value</th></tr></thead><tbody>")
    last_group = None
    for group, metric, value in _summary_rows(report):
        group_cell = e(group) if group != last_group else ""
        cls = " class='group-start'" if group != last_group else ""
        parts.append(f"<tr{cls}><td class='grp'>{group_cell}</td><td>{e(metric)}</td>"
                     f"<td class='val'>{e(str(value))}</td></tr>")
        last_group = group
    parts.append("</tbody></table></section>")

    # Inputs
    if report.inputs:
        parts.append("<section><h2>Inputs</h2><ul class='inputs'>")
        for key, value in report.inputs.items():
            parts.append(f"<li><b>{e(str(key))}:</b> {e(str(value))}</li>")
        parts.append("</ul></section>")

    if report.warnings:
        parts.append("<section><h2>Parse warnings</h2><ul class='warns'>")
        for w in report.warnings:
            parts.append(f"<li>{e(str(w))}</li>")
        parts.append("</ul></section>")

    parts.append(
        "<footer>Thresholds follow the Earth BioGenome Project / VGP "
        "<b>3C</b> framework — Contiguity, Completeness, Correctness. "
        "Generated by AsmQC (GPL-3.0).</footer>"
    )
    parts.append("</div></body></html>")
    return "".join(parts)


def write_html(report: QCReport, path: str | Path, plots: dict[str, str] | None = None) -> Path:
    path = Path(path)
    path.write_text(render_html(report, plots), encoding="utf-8")
    return path


_CSS = """
:root{--fg:#1c1e21;--muted:#6c757d;--line:#e3e6ea;--bg:#f7f8fa;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
 font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;
 line-height:1.5;font-size:15px}
.wrap{max-width:920px;margin:0 auto;padding:24px}
header{background:#fff;border:1px solid var(--line);border-radius:12px;padding:22px;margin-bottom:20px}
h1{margin:0;font-size:24px}
.subtitle{font-size:18px;color:var(--muted);margin-top:2px;word-break:break-all}
.status{display:inline-block;color:#fff;font-weight:600;border-radius:8px;
 padding:6px 14px;margin-top:12px;font-size:15px}
.meta{color:var(--muted);font-size:13px;margin-top:10px}
.chips{margin-top:12px;display:flex;flex-wrap:wrap;gap:6px}
.chip{color:#fff;border-radius:20px;padding:3px 11px;font-size:12px;font-weight:600}
section{background:#fff;border:1px solid var(--line);border-radius:12px;padding:20px;margin-bottom:20px}
h2{margin:0 0 14px;font-size:18px;border-bottom:2px solid var(--line);padding-bottom:8px}
.flag{border:1px solid var(--line);border-radius:10px;padding:14px;margin-bottom:12px;background:#fcfcfd}
.flag-head{display:flex;align-items:center;gap:10px;flex-wrap:wrap}
.sev{color:#fff;font-weight:700;border-radius:6px;padding:2px 9px;font-size:12px;letter-spacing:.3px}
.flag-title{font-weight:600;font-size:16px}
.flag-cat{color:var(--muted);font-size:12px;border:1px solid var(--line);border-radius:12px;padding:1px 9px}
.flag-msg{margin-top:8px}
.flag-why,.flag-act{margin-top:7px;font-size:14px;color:#33363a}
.flag-act{color:#0b5d2e}
.plots{display:grid;grid-template-columns:1fr;gap:18px}
figure{margin:0;text-align:center}
figure img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:8px;background:#fff}
figcaption{color:var(--muted);font-size:13px;margin-top:6px}
table.metrics{width:100%;border-collapse:collapse;font-size:14px}
table.metrics th,table.metrics td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line)}
table.metrics tr.group-start td{border-top:2px solid var(--line)}
td.grp{font-weight:600;color:#0b3a66;white-space:nowrap}
td.val{font-variant-numeric:tabular-nums;text-align:right;white-space:nowrap}
ul.inputs,ul.warns{margin:0;padding-left:18px}
.muted{color:var(--muted)}
footer{color:var(--muted);font-size:13px;text-align:center;padding:8px 0 24px}
@media(min-width:760px){.plots{grid-template-columns:1fr 1fr}}
"""
