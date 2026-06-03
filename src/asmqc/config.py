# AsmQC — configuration loading and deep-merge over built-in defaults.
# Copyright (C) 2026 AsmQC contributors. Licensed under GPL-3.0-or-later.
"""Load curation thresholds.

The built-in defaults live in :file:`data/default_config.yaml`.  A user-supplied
YAML is *deep-merged* over them, so a config file only needs to set the knobs it
wants to change.
"""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

_DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_CONFIG_PATH = _DATA_DIR / "default_config.yaml"


class AsmQCConfig:
    """Thin wrapper around the merged threshold dictionary.

    Access values with dotted helpers::

        cfg.get("busco", "duplicated_warn_pct")     # -> 5.0
        cfg.enabled("contamination")                 # -> True
    """

    def __init__(self, data: dict[str, Any]):
        self.data = data

    # -- construction --------------------------------------------------------
    @classmethod
    def load(cls, path: str | Path | None = None) -> AsmQCConfig:
        """Load defaults, then deep-merge the user config at *path* (if any)."""
        base = _read_yaml(DEFAULT_CONFIG_PATH)
        if path is not None:
            user = _read_yaml(Path(path))
            base = _deep_merge(base, user)
        return cls(base)

    @classmethod
    def default(cls) -> AsmQCConfig:
        return cls(_read_yaml(DEFAULT_CONFIG_PATH))

    # -- access --------------------------------------------------------------
    def section(self, name: str) -> dict[str, Any]:
        value = self.data.get(name, {})
        return value if isinstance(value, dict) else {}

    def get(self, section: str, key: str, default: Any = None) -> Any:
        return self.section(section).get(key, default)

    def enabled(self, family: str) -> bool:
        return bool(self.section("enabled").get(family, True))

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self.data)


def default_config_text() -> str:
    """Return the verbatim default-config YAML (for ``asmqc init-config``)."""
    return DEFAULT_CONFIG_PATH.read_text(encoding="utf-8")


def _read_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        loaded = yaml.safe_load(fh)
    if loaded is None:
        return {}
    if not isinstance(loaded, dict):
        raise ValueError(f"Config root must be a mapping, got {type(loaded).__name__}: {path}")
    return loaded


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge *override* into a copy of *base*."""
    result = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result
