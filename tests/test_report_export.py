# AsmQC tests — report writers, plots and browser-track exporters.
import json

from asmqc.config import AsmQCConfig
from asmqc.core import build_report
from asmqc.exporters import build_features, write_jbrowse_config, write_tracks
from asmqc.models import (
    AssemblyStats,
    BuscoResult,
    ContaminationHit,
    ContaminationResult,
    MerquryResult,
    QCReport,
    TelomereResult,
    TelomereScaffold,
)
from asmqc.plots import build_plots
from asmqc.report import render_html, render_markdown, write_json


def _full_report():
    s = AssemblyStats(source="fasta", total_length=1_000_000, num_sequences=5,
                      num_contigs=6, n50=400_000, n90=100_000, l50=2, contig_n50=200_000,
                      largest=500_000, gc_percent=41.0, gap_count=3, gap_bases=300)
    b = BuscoResult(complete_pct=88, single_pct=80, duplicated_pct=8,
                    fragmented_pct=6, missing_pct=6, lineage="aves_odb10")
    m = MerquryResult(qv=32.0, completeness_pct=83.0)
    t = TelomereResult(repeat_motif="TTAGGG", window_size=10000, chromosome_min_length=100000,
                       scaffolds=[TelomereScaffold("chr1", length=500_000,
                                                   start_telomere=True, end_telomere=True)])
    c = ContaminationResult(source="fcs", n_sequences_flagged=1, total_flagged_bases=76092,
                            hits=[ContaminationHit("contig_017", action="EXCLUDE",
                                                   start=1, end=76092, taxon="Ponticoccus",
                                                   note="prok:a-proteobacteria")])
    r = QCReport("test_asm", stats=s, busco=b, merqury=m, telomere=t, contamination=c,
                 asmqc_version="0.1.0", generated_at="2026-06-03 00:00:00 UTC")
    from asmqc.flags import evaluate
    r.flags = evaluate(r, AsmQCConfig.default())
    return r


def test_json_serialisable(tmp_path):
    r = _full_report()
    p = write_json(r, tmp_path / "summary.json")
    data = json.loads(p.read_text())
    assert data["assembly_name"] == "test_asm"
    assert data["overall_status"] == "FLAG"
    assert data["busco"]["duplicated_pct"] == 8.0
    assert data["assembly_stats"]["n50"] == 400000
    assert isinstance(data["flags"], list) and data["flags"]


def test_markdown_renders():
    md = render_markdown(_full_report())
    assert "# AsmQC report" in md
    assert "Curation flags" in md
    assert "| Group | Metric | Value |" in md
    assert "Scaffold N50" in md


def test_html_self_contained():
    r = _full_report()
    plots = build_plots(r, AsmQCConfig.default())
    html = render_html(r, plots)
    assert html.startswith("<!DOCTYPE html>")
    assert "AsmQC report" in html
    # Plots embedded as base64 data URIs (self-contained, no external assets).
    assert "data:image/png;base64," in html
    assert "http://" not in html  # no remote asset references


def test_plots_built():
    r = _full_report()
    plots = build_plots(r, AsmQCConfig.default())
    assert "busco" in plots
    # Each plot value is non-empty base64.
    assert all(isinstance(v, str) and len(v) > 100 for v in plots.values())


def test_tracks_export(tmp_path, fixtures):
    cfg = AsmQCConfig.default()
    r = build_report(fasta=fixtures / "tiny.fasta",
                     contamination=fixtures / "fcs_gx_report.txt",
                     tidk=fixtures / "tidk_windows.tsv", config=cfg)
    features = build_features(r, fasta_path=fixtures / "tiny.fasta", min_gap_len=10)
    ftypes = {f.ftype for f in features}
    assert "gap" in ftypes          # from FASTA
    assert "telomere" in ftypes     # from tidk
    assert "contamination" in ftypes  # from FCS

    out = write_tracks(r, tmp_path, fasta_path=fixtures / "tiny.fasta")
    gff = out["gff3"].read_text()
    assert gff.startswith("##gff-version 3")
    assert "\tgap\t" in gff
    bed = out["bed"].read_text()
    assert "itemRgb" in bed


def test_gff_bed_coordinates(tmp_path):
    # BED start must be GFF start - 1 (0-based vs 1-based).
    from asmqc.exporters.tracks import TrackFeature, write_bed, write_gff3
    feat = [TrackFeature("chr1", 101, 200, "gap", "gap_1")]
    g = write_gff3(feat, tmp_path / "f.gff3").read_text().strip().splitlines()[-1].split("\t")
    b = write_bed(feat, tmp_path / "f.bed").read_text().strip().splitlines()[-1].split("\t")
    assert g[3] == "101" and g[4] == "200"
    assert b[1] == "100" and b[2] == "200"


def test_gff_bed_seqid_escaping(tmp_path):
    # A space/tab in a seq_id must not shift the tab-delimited columns.
    from asmqc.exporters.tracks import TrackFeature, write_bed, write_gff3
    feat = [TrackFeature("scaffold 1|arrow", 1, 100, "contamination", "EXCLUDE_1")]
    gff = write_gff3(feat, tmp_path / "f.gff3").read_text().strip().splitlines()[-1]
    bed = write_bed(feat, tmp_path / "f.bed").read_text().strip().splitlines()[-1]
    # GFF3 row still has exactly 9 tab columns and the seqid has no raw space.
    gcols = gff.split("\t")
    assert len(gcols) == 9
    assert " " not in gcols[0] and gcols[0] != "scaffold 1|arrow"
    bcols = bed.split("\t")
    assert len(bcols) == 9
    assert " " not in bcols[0]


def test_markdown_escapes_table_breaking_chars():
    # A pipe / markup in an untrusted field must not corrupt the table or inject HTML.
    from asmqc.models import BuscoResult
    r = QCReport("asm|<b>", busco=BuscoResult(complete_pct=99, lineage="aves|odb<script>"))
    md = render_markdown(r)
    # The injected pipe is escaped, and angle brackets are neutralised.
    assert "aves\\|odb" in md
    assert "<script>" not in md


def test_jbrowse_config(tmp_path):
    p = write_jbrowse_config(tmp_path, "test_asm", "asm.fa", gff3_name="flags.gff3")
    cfg = json.loads(p.read_text())
    assert cfg["assemblies"][0]["name"] == "test_asm"
    adapter = cfg["assemblies"][0]["sequence"]["adapter"]
    assert adapter["type"] == "IndexedFastaAdapter"
    assert adapter["fastaLocation"]["uri"] == "asm.fa"
    assert cfg["tracks"][0]["adapter"]["type"] == "Gff3Adapter"
