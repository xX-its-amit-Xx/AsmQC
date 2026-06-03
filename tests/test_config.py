# AsmQC tests — configuration loading and merge.
import pytest

from asmqc.config import AsmQCConfig, default_config_text


def test_default_loads():
    cfg = AsmQCConfig.default()
    assert cfg.get("merqury", "qv_warn") == 40.0
    assert cfg.get("busco", "duplicated_warn_pct") == 5.0
    assert cfg.enabled("contamination") is True


def test_deep_merge_override(tmp_path):
    p = tmp_path / "ov.yaml"
    p.write_text("merqury:\n  qv_warn: 50\nenabled:\n  busco: false\n")
    cfg = AsmQCConfig.load(p)
    # Overridden value.
    assert cfg.get("merqury", "qv_warn") == 50
    # Sibling default preserved (not wiped by the partial override).
    assert cfg.get("merqury", "qv_flag") == 30.0
    assert cfg.enabled("busco") is False
    # Untouched section intact.
    assert cfg.get("busco", "duplicated_warn_pct") == 5.0


def test_init_config_text_roundtrips(tmp_path):
    text = default_config_text()
    p = tmp_path / "c.yaml"
    p.write_text(text)
    cfg = AsmQCConfig.load(p)
    assert cfg.get("contiguity", "scaffold_n50_warn_bp") == 1000000


def test_missing_config_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        AsmQCConfig.load(tmp_path / "missing.yaml")
