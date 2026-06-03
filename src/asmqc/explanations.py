# AsmQC — plain-English explanations and curation actions for each flag.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Curator-facing prose attached to every flag.

Each entry maps a logical issue key to a likely *cause* (``explanation``) and a
recommended *curation action* (``action``).  The flag engine fills the numeric
``message`` separately; this module is the "why / what next" knowledge base so a
non-expert can act on a flag without reading a methods paper.
"""
from __future__ import annotations

EXPLANATIONS: dict[str, dict[str, str]] = {
    # -- BUSCO ---------------------------------------------------------------
    "busco_duplication": {
        "explanation": (
            "BUSCO looks for genes that should occur exactly once in this lineage. "
            "When many of them appear twice, the assembly is probably carrying both "
            "haplotypes of a diploid genome in the same (supposedly haploid) set of "
            "sequences — i.e. uncollapsed 'haplotigs' / false duplications. A truly "
            "duplicated region (recent whole-genome duplication) is possible but rarer."
        ),
        "action": (
            "Run purge_dups (or purge_haplotigs) to collapse retained haplotigs into a "
            "primary haplotype, then re-check BUSCO duplication and the Merqury copy-number "
            "spectrum (the read-only/2-copy peaks). If duplication is biological, document "
            "the expected ploidy/WGD."
        ),
    },
    "busco_completeness": {
        "explanation": (
            "A low Complete percentage means expected single-copy genes are not found "
            "intact in the assembly. Causes: genuinely missing sequence (low coverage, "
            "collapsed repeats), genes broken across contig boundaries, or — for distant "
            "lineages — real biological gene loss. High Missing+Fragmented usually points "
            "to a fragmented or incomplete assembly rather than biology."
        ),
        "action": (
            "Check coverage and k-mer completeness (Merqury) to see if sequence is missing; "
            "consider additional long-read data or a re-assembly. Confirm the BUSCO lineage "
            "dataset is appropriate for the species before concluding genes are truly absent."
        ),
    },
    "busco_fragmented": {
        "explanation": (
            "Fragmented BUSCOs are found only partially — typically a gene split across two "
            "contigs/scaffolds because the assembly breaks inside it. This tracks low "
            "contiguity: short contigs interrupt gene models."
        ),
        "action": (
            "Improve contiguity (long-read scaffolding, gap filling). Fragmentation usually "
            "falls as contig N50 rises; re-evaluate after scaffolding."
        ),
    },
    "busco_missing": {
        "explanation": (
            "Missing BUSCOs are not detected at all. Combined with low completeness this "
            "indicates absent sequence; alone it can also reflect an over-stringent or "
            "mismatched lineage dataset."
        ),
        "action": (
            "Verify the lineage dataset, then investigate coverage gaps. Persistent high "
            "Missing after good coverage suggests real dropout regions to target with "
            "additional data."
        ),
    },
    # -- Merqury -------------------------------------------------------------
    "merqury_qv": {
        "explanation": (
            "Merqury's QV is a reference-free estimate of per-base consensus accuracy from "
            "k-mers that appear in the assembly but not the reads (assumed errors). A low QV "
            "means residual base-level errors — common with noisy long reads that were not "
            "polished, or under-polished consensus. Q40 ≈ 1 error / 10 kb (the EBP standard); "
            "below Q30 most genes will carry a base error."
        ),
        "action": (
            "Polish the consensus (e.g. with the HiFi reads via DeepVariant/racon/medaka as "
            "appropriate) and re-run Merqury. Check the k-mer spectrum for a residual error "
            "peak at low multiplicity."
        ),
    },
    "merqury_completeness": {
        "explanation": (
            "K-mer completeness is the fraction of reliable read k-mers that are present in "
            "the assembly. Below ~90% means real genomic sequence (often repeats or one "
            "haplotype) is missing from the assembly even if gene space looks fine."
        ),
        "action": (
            "Investigate missing-sequence sources: collapsed repeats, dropped haplotype, or "
            "low coverage regions. Additional long reads or a re-assembly may recover the "
            "missing k-mers."
        ),
    },
    # -- Contiguity ----------------------------------------------------------
    "contiguity_scaffold_n50": {
        "explanation": (
            "Scaffold N50 summarises how long the scaffolds are: half the genome sits in "
            "scaffolds at least this long. A low value means the assembly has not reached "
            "chromosome scale — it is a fragmented draft rather than a reference. The EBP "
            "target is chromosome-scale scaffolds (often 10 Mb+ N50)."
        ),
        "action": (
            "Add long-range information (Hi-C, optical maps, or linkage) to scaffold contigs "
            "into chromosomes, then re-assess. Confirm against the expected karyotype."
        ),
    },
    "contiguity_contig_n50": {
        "explanation": (
            "Contig N50 measures the underlying contiguity before scaffolding (gaps removed). "
            "Short contigs mean the primary assembly is broken — usually short or noisy reads, "
            "or unresolved repeats. Scaffolding can join these but cannot recover the missing "
            "sequence in the gaps."
        ),
        "action": (
            "Prefer long high-accuracy reads (PacBio HiFi / ONT) and re-assemble; contig N50 "
            "is set at the contigging stage, not by scaffolding."
        ),
    },
    "contiguity_fragmentation": {
        "explanation": (
            "A large share of the assembly sits in many short sequences. This 'debris' is "
            "typically unplaced fragments, alternate haplotigs, or contamination, and it "
            "inflates the sequence count while contributing little to the chromosomes."
        ),
        "action": (
            "Purge duplicates, screen for contamination (FCS/Kraken), and consider removing or "
            "binning very short unplaced scaffolds. Investigate why so much sequence stayed "
            "unscaffolded."
        ),
    },
    "contiguity_many_sequences": {
        "explanation": (
            "The assembly contains a very large number of sequences for its size. On its own "
            "this is a smell rather than a verdict, but it usually accompanies fragmentation, "
            "retained haplotigs, or contaminant contigs."
        ),
        "action": (
            "Cross-check with the fragmentation, duplication and contamination flags; reduce "
            "the unplaced/short fraction during curation."
        ),
    },
    # -- Gaps ----------------------------------------------------------------
    "gaps_density": {
        "explanation": (
            "Gaps are runs of N inserted where scaffolding joined contigs without knowing the "
            "intervening sequence. Many gaps mean the chromosome-scale structure rests on "
            "unfilled joins; each gap is also a place a structural error can hide. A "
            "telomere-to-telomere assembly has zero gaps."
        ),
        "action": (
            "Gap-fill with long reads (e.g. ONT/HiFi targeted assembly) and validate joins "
            "with Hi-C. Track gaps-per-Gbp toward the EBP target (≤1000/Gbp)."
        ),
    },
    # -- Telomeres -----------------------------------------------------------
    "telomere_incomplete": {
        "explanation": (
            "Telomeric repeat arrays were not detected at both ends of most chromosome-scale "
            "scaffolds. That means the assembly is not yet telomere-to-telomere: some "
            "chromosome ends are missing or unresolved (telomeres are highly repetitive and "
            "hard to assemble). This is expected for most non-T2T assemblies and is a 'note', "
            "not a failure."
        ),
        "action": (
            "If T2T completeness is a goal, target the missing ends with ultra-long reads and "
            "verify with the telomere map. Otherwise record which chromosomes are arm-to-arm "
            "complete."
        ),
    },
    "telomere_low": {
        "explanation": (
            "Telomeric repeats were found at very few chromosome ends. Either the canonical "
            "motif searched does not match this clade's telomere repeat, or chromosome ends "
            "are largely unassembled."
        ),
        "action": (
            "Confirm the correct telomere motif for the clade (tidk explore / find) before "
            "concluding ends are missing; then target unresolved ends with long reads."
        ),
    },
    # -- Contamination -------------------------------------------------------
    "contamination_fcs": {
        "explanation": (
            "NCBI FCS-GX flagged sequence whose best taxonomic match is a foreign organism "
            "(bacteria, other host, adaptor, organelle in the wrong place). EXCLUDE/TRIM "
            "actions are high-confidence foreign sequence that NCBI would remove before "
            "release; it often hides in small, low-coverage contigs."
        ),
        "action": (
            "Remove EXCLUDE sequences and TRIM the flagged spans (the GFF/BED export marks "
            "exact coordinates). Re-run FCS after cleaning and confirm no on-target sequence "
            "was lost. REVIEW calls need a manual look."
        ),
    },
    "contamination_kraken": {
        "explanation": (
            "Kraken2 classified a non-trivial fraction of the assembly's k-mers to a "
            "superkingdom other than the target organism's — a sign of contaminant sequence "
            "(e.g. bacterial contigs in a eukaryotic assembly) or, less often, mislabelled "
            "reads. Kraken is read/k-mer level, so it flags presence, not exact coordinates."
        ),
        "action": (
            "Confirm with FCS-GX or BLAST for per-sequence coordinates, then remove "
            "contaminant contigs. Check whether the contaminant correlates with the short "
            "unplaced scaffold fraction."
        ),
    },
    # -- All-clear -----------------------------------------------------------
    "pass": {
        "explanation": "This metric meets the configured quality standard.",
        "action": "No action required.",
    },
}


def get(key: str) -> tuple[str, str]:
    """Return ``(explanation, action)`` for an issue key."""
    entry = EXPLANATIONS.get(key, EXPLANATIONS["pass"])
    return entry["explanation"], entry["action"]
