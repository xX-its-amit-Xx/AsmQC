# AsmQC tests — tool output parsers.
import pytest

from asmqc.parsers import (
    parse_busco,
    parse_contamination,
    parse_gfastats,
    parse_merqury,
    parse_quast,
    parse_tidk,
)
from asmqc.parsers.contamination import parse_fcs, parse_kraken


# --------------------------------------------------------------------------- #
# BUSCO
# --------------------------------------------------------------------------- #
def test_busco_txt(fixtures):
    b = parse_busco(fixtures / "busco_short_summary.txt")
    assert b is not None
    assert b.complete_pct == 88.0
    assert b.single_pct == 80.0
    assert b.duplicated_pct == 8.0
    assert b.fragmented_pct == 6.0
    assert b.missing_pct == 6.0
    assert b.complete_n == 7338
    assert b.duplicated_n == 667
    assert b.total_n == 8338
    assert b.lineage == "aves_odb10"
    assert b.busco_version == "5.4.3"
    assert b.mode == "genome"


def test_busco_json(fixtures):
    b = parse_busco(fixtures / "busco_short_summary.json")
    assert b is not None
    assert b.complete_pct == 97.6
    assert b.duplicated_pct == 2.6  # "Multi copy percentage"
    assert b.single_n == 13091
    assert b.duplicated_n == 358
    assert b.total_n == 13780
    assert b.lineage == "primates_odb10"
    assert b.busco_version == "6.0.0"


def test_busco_counts_only(tmp_path):
    p = tmp_path / "s.txt"
    p.write_text(
        "\t100\tComplete BUSCOs (C)\n"
        "\t90\tComplete and single-copy BUSCOs (S)\n"
        "\t10\tComplete and duplicated BUSCOs (D)\n"
        "\t3\tFragmented BUSCOs (F)\n"
        "\t2\tMissing BUSCOs (M)\n"
        "\t105\tTotal BUSCO groups searched\n"
    )
    b = parse_busco(p)
    assert b.complete_n == 100
    # Percentages derived from counts.
    assert b.complete_pct == pytest.approx(100 / 105 * 100, abs=0.01)


def test_busco_missing_file(tmp_path):
    assert parse_busco(tmp_path / "nope.txt") is None


# --------------------------------------------------------------------------- #
# Merqury
# --------------------------------------------------------------------------- #
def test_merqury_dir(fixtures):
    m = parse_merqury(fixtures / "merqury")
    assert m is not None
    assert m.qv == pytest.approx(35.1183)
    assert m.kmers_asm_only == 682359
    assert m.kmers_total == 123511626
    assert m.completeness_pct == pytest.approx(83.7764)
    assert len(m.spectra) == 10
    classes = {p.copy_class for p in m.spectra}
    assert "read-only" in classes and ">4" in classes


def test_merqury_qv_only(fixtures):
    m = parse_merqury(qv=fixtures / "merqury" / "asm.qv")
    assert m.qv == pytest.approx(35.1183)
    assert m.completeness_pct is None


def test_merqury_both_row(tmp_path):
    p = tmp_path / "x.qv"
    p.write_text("hap1\t10\t1000\t30.0\t0.001\nhap2\t12\t1000\t29.0\t0.001\nBoth\t20\t2000\t31.0\t0.0008\n")
    m = parse_merqury(qv=p)
    assert m.qv == 31.0  # 'Both' row preferred
    assert set(m.per_assembly_qv) == {"hap1", "hap2", "Both"}


# --------------------------------------------------------------------------- #
# QUAST / gfastats
# --------------------------------------------------------------------------- #
def test_quast(fixtures):
    s = parse_quast(fixtures / "quast_report.tsv")
    assert s is not None
    assert s.source == "quast"
    assert s.total_length == 1026755127
    assert s.n50 == 71229700
    assert s.n90 == 22000000
    assert s.l50 == 5
    assert s.gc_percent == 41.5
    assert s.num_sequences == 164
    assert s.n_per_100kbp == 7.01


def test_gfastats(fixtures):
    s = parse_gfastats(fixtures / "gfastats_summary.txt")
    assert s is not None
    assert s.source == "gfastats"
    assert s.num_scaffolds == 164
    assert s.total_length == 1026755127
    assert s.n50 == 71229700
    assert s.l50 == 5
    assert s.num_contigs == 374
    assert s.contig_n50 == 17905263
    assert s.gap_count == 210
    assert s.gap_bases == 71937
    assert s.gc_percent == 41.5


# --------------------------------------------------------------------------- #
# tidk
# --------------------------------------------------------------------------- #
def test_tidk(fixtures):
    t = parse_tidk(fixtures / "tidk_windows.tsv", min_window_count=50, terminal_windows=2)
    assert t is not None
    assert t.window_size == 10000
    assert t.repeat_motif == "TTAGGG"
    assert t.n_scaffolds == 3
    assert t.n_with_both == 1   # chr1
    assert t.n_with_any == 2    # chr1, chr2
    chr1 = next(s for s in t.scaffolds if s.name == "chr1")
    assert chr1.both_ends is True
    chr3 = next(s for s in t.scaffolds if s.name == "chr3")
    assert chr3.any_end is False


def test_tidk_motif_not_poisoned_by_leading_blank(tmp_path):
    # A leading blank line pushes the header off line 0; the motif must still be
    # the real repeat, not the 'telomeric_repeat' column header.
    p = tmp_path / "w.tsv"
    p.write_text(
        "\nid\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\n"
        "chr1\t10000\t3\t412\tTTAGGG\nchr1\t50000\t389\t2\tTTAGGG\n"
    )
    t = parse_tidk(p, min_window_count=50, terminal_windows=2)
    assert t.repeat_motif == "TTAGGG"


# --------------------------------------------------------------------------- #
# Contamination
# --------------------------------------------------------------------------- #
def test_fcs(fixtures):
    c = parse_fcs(fixtures / "fcs_gx_report.txt")
    assert c is not None
    assert c.source == "fcs"
    assert c.n_hits == 3
    actions = {h.action for h in c.hits}
    assert actions == {"EXCLUDE", "TRIM", "REVIEW"}
    assert c.n_sequences_flagged == 3
    # EXCLUDE span 1..76092 = 76092, TRIM 500..2300 = 1801, REVIEW 1..1200 = 1200.
    assert c.total_flagged_bases == 76092 + 1801 + 1200


def test_kraken(fixtures):
    c = parse_kraken(fixtures / "kraken_report.txt")
    assert c is not None
    assert c.source == "kraken"
    assert c.unclassified_pct == 10.0
    assert c.top_taxa.get("Eukaryota") == 85.0
    assert c.top_taxa.get("Bacteria") == 4.5


def test_contamination_autodetect(fixtures):
    assert parse_contamination(fixtures / "fcs_gx_report.txt").source == "fcs"
    assert parse_contamination(fixtures / "kraken_report.txt").source == "kraken"
