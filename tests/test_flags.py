# AsmQC tests — curation flag engine.
from asmqc.config import AsmQCConfig
from asmqc.flags import evaluate
from asmqc.models import (
    AssemblyStats,
    BuscoResult,
    ContaminationHit,
    ContaminationResult,
    MerquryResult,
    QCReport,
    SeqStat,
    Severity,
    TelomereResult,
    TelomereScaffold,
)


def _cfg():
    return AsmQCConfig.default()


def _flags_by_id(report):
    return {f.id: f for f in report.flags}


def test_busco_duplication_warn_and_flag():
    cfg = _cfg()
    r = QCReport("a", busco=BuscoResult(complete_pct=96, duplicated_pct=8, single_pct=88))
    r.flags = evaluate(r, cfg)
    f = _flags_by_id(r)["busco_duplication"]
    assert f.severity == Severity.WARN
    assert "haplotype" in f.explanation.lower()

    r2 = QCReport("a", busco=BuscoResult(complete_pct=96, duplicated_pct=12, single_pct=84))
    r2.flags = evaluate(r2, cfg)
    assert _flags_by_id(r2)["busco_duplication"].severity == Severity.FLAG


def test_busco_completeness_flag():
    r = QCReport("a", busco=BuscoResult(complete_pct=85, duplicated_pct=1))
    r.flags = evaluate(r, _cfg())
    assert _flags_by_id(r)["busco_completeness"].severity == Severity.FLAG


def test_busco_pass():
    r = QCReport("a", busco=BuscoResult(complete_pct=98, duplicated_pct=1,
                                        fragmented_pct=0.5, missing_pct=0.5))
    r.flags = evaluate(r, _cfg())
    assert _flags_by_id(r)["busco_ok"].severity == Severity.PASS


def test_merqury_qv_thresholds():
    r = QCReport("a", merqury=MerquryResult(qv=32))
    r.flags = evaluate(r, _cfg())
    assert _flags_by_id(r)["merqury_qv"].severity == Severity.WARN

    r2 = QCReport("a", merqury=MerquryResult(qv=25))
    r2.flags = evaluate(r2, _cfg())
    assert _flags_by_id(r2)["merqury_qv"].severity == Severity.FLAG

    r3 = QCReport("a", merqury=MerquryResult(qv=45))
    r3.flags = evaluate(r3, _cfg())
    assert _flags_by_id(r3)["merqury_qv"].severity == Severity.PASS


def test_contiguity_scaffold_n50():
    s = AssemblyStats(source="fasta", total_length=10_000_000, num_sequences=10, n50=500_000)
    r = QCReport("a", stats=s)
    r.flags = evaluate(r, _cfg())
    assert _flags_by_id(r)["contiguity_scaffold_n50"].severity == Severity.WARN


def test_fragmentation_flag():
    per_seq = [SeqStat(name="big", length=900_000)]
    per_seq += [SeqStat(name=f"s{i}", length=5000) for i in range(40)]  # 200k in short
    total = 900_000 + 40 * 5000
    s = AssemblyStats(source="fasta", total_length=total, num_sequences=41,
                      n50=900_000, per_seq=per_seq, num_short_seqs=40)
    r = QCReport("a", stats=s)
    r.flags = evaluate(r, _cfg())
    assert "contiguity_fragmentation" in _flags_by_id(r)


def test_gap_density_note():
    s = AssemblyStats(source="fasta", total_length=1_000_000, num_sequences=1,
                      gap_count=100, gap_bases=2000)
    r = QCReport("a", stats=s)
    r.flags = evaluate(r, _cfg())
    # 100 gaps / 1 Mb = 10 per 100 kb >= 5 -> note (gap fraction 0.2% < 5% warn).
    assert _flags_by_id(r)["gaps_density"].severity == Severity.NOTE


def test_telomere_not_t2t():
    scaffolds = [
        TelomereScaffold("chr1", length=2_000_000, start_telomere=True, end_telomere=True),
        TelomereScaffold("chr2", length=2_000_000, start_telomere=True, end_telomere=False),
    ]
    t = TelomereResult(repeat_motif="TTAGGG", window_size=10000, scaffolds=scaffolds,
                       chromosome_min_length=1_000_000)
    r = QCReport("a", telomere=t)
    r.flags = evaluate(r, _cfg())
    assert _flags_by_id(r)["telomere_t2t"].severity == Severity.NOTE


def test_contamination_fcs_flag():
    c = ContaminationResult(source="fcs", n_sequences_flagged=1, total_flagged_bases=76092,
                            hits=[ContaminationHit("c1", action="EXCLUDE", start=1, end=76092,
                                                   taxon="Ponticoccus")])
    r = QCReport("a", contamination=c)
    r.flags = evaluate(r, _cfg())
    assert _flags_by_id(r)["contamination_fcs"].severity == Severity.FLAG


def test_contamination_fcs_totals_exclude_review():
    # The FLAG message/evidence must total only the EXCLUDE/TRIM subset, not the
    # parser's all-action totals (which include REVIEW/INFO).
    c = ContaminationResult(
        source="fcs",
        n_sequences_flagged=3, total_flagged_bases=79093,  # parser-wide (incl. REVIEW)
        hits=[
            ContaminationHit("c1", action="EXCLUDE", start=1, end=76092),
            ContaminationHit("c2", action="TRIM", start=500, end=2300),
            ContaminationHit("c3", action="REVIEW", start=1, end=1200),
        ],
    )
    r = QCReport("a", contamination=c)
    r.flags = evaluate(r, _cfg())
    fcs = _flags_by_id(r)["contamination_fcs"]
    # EXCLUDE 76092 + TRIM 1801 = 77893, across 2 sequences (REVIEW excluded).
    assert fcs.evidence["total_flagged_bases"] == 77893
    assert fcs.evidence["n_sequences_flagged"] == 2
    assert "77,893 bp" in fcs.message


def test_contamination_taxon_sanitised():
    c = ContaminationResult(
        source="fcs", n_sequences_flagged=1, total_flagged_bases=100,
        hits=[ContaminationHit("c1", action="EXCLUDE", start=1, end=100,
                               taxon="<script>x</script>|evil")],
    )
    r = QCReport("a", contamination=c)
    r.flags = evaluate(r, _cfg())
    msg = _flags_by_id(r)["contamination_fcs"].message
    assert "<script>" not in msg and "|evil" not in msg


def test_contamination_kraken_flag():
    c = ContaminationResult(source="kraken", top_taxa={"Eukaryota": 85.0, "Bacteria": 4.5},
                            unclassified_pct=10.0)
    r = QCReport("a", contamination=c)
    r.flags = evaluate(r, _cfg())
    assert _flags_by_id(r)["contamination_kraken"].severity == Severity.FLAG


def test_disabled_family_silences_flags():
    cfg = AsmQCConfig.load(None)
    cfg.data["enabled"]["busco"] = False
    r = QCReport("a", busco=BuscoResult(complete_pct=50, duplicated_pct=30))
    r.flags = evaluate(r, cfg)
    assert all(not f.id.startswith("busco") for f in r.flags)


def test_overall_severity_is_max():
    r = QCReport("a", busco=BuscoResult(complete_pct=85, duplicated_pct=1),
                 merqury=MerquryResult(qv=45))
    r.flags = evaluate(r, _cfg())
    assert r.overall_severity == Severity.FLAG
