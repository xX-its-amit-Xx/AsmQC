# AsmQC tests — FASTA contiguity stats.
from asmqc.fasta_stats import compute_assembly_stats, iter_fasta, nx_lx


def test_iter_fasta_ids_and_lengths(fixtures):
    records = dict(iter_fasta(fixtures / "tiny.fasta"))
    # Header description is stripped to the first token.
    assert set(records) == {"seqA", "seqB", "seqC", "seqD", "seqE"}
    assert len(records["seqA"]) == 100
    assert len(records["seqB"]) == 80
    assert len(records["seqC"]) == 60
    assert len(records["seqE"]) == 20


def test_basic_contiguity(fixtures):
    s = compute_assembly_stats(fixtures / "tiny.fasta", min_gap_len=10)
    assert s.total_length == 300
    assert s.num_sequences == 5
    assert s.largest == 100
    assert s.smallest == 20
    # Sorted desc 100,80,60,40,20: N50 reaches 150 at 80 (idx 2).
    assert s.n50 == 80
    assert s.l50 == 2
    # N90: 270 reached at the 40 bp seq (idx 4).
    assert s.n90 == 40
    assert s.l90 == 4


def test_gaps_and_contigs(fixtures):
    s = compute_assembly_stats(fixtures / "tiny.fasta", min_gap_len=10)
    # seqC has one 10-N gap -> 1 gap, 10 gap bases, splits into 2 contigs.
    assert s.gap_count == 1
    assert s.gap_bases == 10
    # Contigs: A100, B80, C->25+25, D40, E20 = 6 contigs.
    assert s.num_contigs == 6
    # Contig lengths desc: 100,80,40,25,25,20 (total 290); half=145 -> contig N50 = 80.
    assert s.contig_n50 == 80
    # N per 100 kbp = 10 / 300 * 1e5.
    assert abs(s.n_per_100kbp - (10 / 300 * 1e5)) < 1e-6


def test_gc_excludes_n(fixtures):
    s = compute_assembly_stats(fixtures / "tiny.fasta", min_gap_len=10)
    # GC = 120 (80 C + 40 G), AT = 170, denom 290.
    assert abs(s.gc_percent - (120 / 290 * 100)) < 1e-6


def test_short_seq_counting(fixtures):
    s = compute_assembly_stats(fixtures / "tiny.fasta", short_scaffold_length=50)
    # seqD (40) and seqE (20) are < 50.
    assert s.num_short_seqs == 2


def test_nx_lx_helper():
    lengths = [100, 80, 60, 40, 20]
    assert nx_lx(lengths, 300, 0.5) == (80, 2)
    assert nx_lx([], 0, 0.5) == (None, None)


def test_empty_fasta(tmp_path):
    p = tmp_path / "empty.fa"
    p.write_text("")
    s = compute_assembly_stats(p)
    assert s.total_length == 0
    assert s.num_sequences == 0
    assert s.n50 is None


def test_pure_n_scaffold_is_gap_not_contig(tmp_path):
    # A whole-record N run must count only as a gap, never also as a contig
    # (else its N bases are double-counted in gap_bases AND contig_n50).
    p = tmp_path / "spacer.fa"
    p.write_text(">realA\n" + "ACGTACGTAC\n" + ">realB\n" + "ACGTACGTAC\n"
                 + ">spacer\n" + "N" * 5000 + "\n")
    s = compute_assembly_stats(p, min_gap_len=10)
    # Two real 10 bp contigs only; the spacer is a single gap.
    assert s.num_contigs == 2
    assert s.contig_n50 == 10
    assert s.gap_count == 1
    assert s.gap_bases == 5000
    spacer = next(x for x in s.per_seq if x.name == "spacer")
    assert spacer.num_subcontigs == 0


def test_gzip_roundtrip(tmp_path):
    import gzip
    p = tmp_path / "a.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">x\n" + "ACGT" * 25 + "\n")
    s = compute_assembly_stats(p)
    assert s.total_length == 100
    assert s.num_sequences == 1
