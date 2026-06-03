# AsmQC tests — shared fixtures.
# Licensed under GPL-3.0-or-later.
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures() -> Path:
    return FIXTURES
