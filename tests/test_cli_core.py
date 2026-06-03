# AsmQC tests — CLI and core orchestration (incl. graceful degradation).
import json

from asmqc.cli import main
from asmqc.config import AsmQCConfig
from asmqc.core import build_report


def test_build_report_partial_inputs(fixtures):
    # Only a FASTA -> stats + contiguity flags, nothing else, no crash.
    r = build_report(fasta=fixtures / "tiny.fasta", config=AsmQCConfig.default())
    assert r.stats is not None
    assert r.busco is None and r.merqury is None
    assert r.telomere is None and r.contamination is None
    # Flags still computed from contiguity alone.
    assert isinstance(r.flags, list)


def test_build_report_no_inputs():
    r = build_report(config=AsmQCConfig.default(), assembly_name="empty")
    assert r.stats is None
    assert r.flags == []
    assert r.overall_severity.name == "PASS"


def test_build_report_all_inputs(fixtures):
    r = build_report(
        fasta=fixtures / "tiny.fasta",
        busco=fixtures / "busco_short_summary.txt",
        merqury=fixtures / "merqury",
        tidk=fixtures / "tidk_windows.tsv",
        contamination=fixtures / "fcs_gx_report.txt",
        config=AsmQCConfig.default(),
    )
    assert r.busco and r.merqury and r.telomere and r.contamination
    ids = {f.id for f in r.flags}
    assert "busco_duplication" in ids
    assert "contamination_fcs" in ids


def test_quast_overrides_fasta_n50(fixtures):
    # When both FASTA and QUAST are supplied, headline N50 comes from QUAST.
    r = build_report(fasta=fixtures / "tiny.fasta", quast=fixtures / "quast_report.tsv",
                     config=AsmQCConfig.default())
    assert r.stats.n50 == 71229700
    assert r.stats.source.startswith("fasta+quast")
    # per-seq richness from FASTA retained.
    assert r.stats.per_seq


def test_malformed_input_warns(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("this is not a busco summary at all\n\x00\x01")
    r = build_report(busco=bad, config=AsmQCConfig.default())
    # No crash; busco simply unparsed.
    assert r.busco is None


def test_cli_run_writes_outputs(fixtures, tmp_path):
    out = tmp_path / "report"
    code = main([
        "run",
        "--fasta", str(fixtures / "tiny.fasta"),
        "--busco", str(fixtures / "busco_short_summary.txt"),
        "--merqury", str(fixtures / "merqury"),
        "--tidk", str(fixtures / "tidk_windows.tsv"),
        "--contamination", str(fixtures / "fcs_gx_report.txt"),
        "--out", str(out),
        "--tracks", "--jbrowse",
        "--quiet",
    ])
    assert code == 0
    assert (out / "report.html").exists()
    assert (out / "report.md").exists()
    assert (out / "summary.json").exists()
    assert (out / "flags.gff3").exists()
    assert (out / "flags.bed").exists()
    assert (out / "jbrowse2_config.json").exists()
    data = json.loads((out / "summary.json").read_text())
    assert data["overall_status"] == "FLAG"


def test_cli_fail_on_gates_exit_code(fixtures, tmp_path):
    base = ["run", "--fasta", str(fixtures / "tiny.fasta"),
            "--contamination", str(fixtures / "fcs_gx_report.txt"),
            "--out", str(tmp_path / "r"), "--quiet"]
    # FCS EXCLUDE -> overall FLAG. --fail-on FLAG should exit 1.
    assert main([*base, "--fail-on", "FLAG"]) == 1
    # A bar above the actual status stays 0.
    assert main([*base, "--fail-on", "FAIL"]) == 0
    # Without --fail-on, a flagged report is still a successful run (exit 0).
    assert main(base) == 0


def test_cli_requires_input(tmp_path):
    code = main(["run", "--out", str(tmp_path / "r")])
    assert code == 2


def test_cli_init_config(capsys):
    code = main(["init-config"])
    assert code == 0
    out = capsys.readouterr().out
    assert "qv_warn" in out
    assert "AsmQC default configuration" in out


def test_cli_version(capsys):
    code = main(["version"])
    assert code == 0
    assert "AsmQC" in capsys.readouterr().out
