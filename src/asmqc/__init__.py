# AsmQC — aggregate genome-assembly QC outputs into one interpretable report.
# Copyright (C) 2026 AsmQC contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""AsmQC: genome-assembly quality-control report aggregator.

Public API::

    from asmqc import run_report, AsmQCConfig
    report = run_report(fasta="asm.fa", busco="short_summary.txt", out_dir="report/")
"""

__version__ = "0.1.0"

from asmqc.models import (  # noqa: E402,F401
    AssemblyStats,
    BuscoResult,
    ContaminationResult,
    Flag,
    MerquryResult,
    QCReport,
    Severity,
    TelomereResult,
)

__all__ = [
    "__version__",
    "AssemblyStats",
    "BuscoResult",
    "ContaminationResult",
    "Flag",
    "MerquryResult",
    "QCReport",
    "Severity",
    "TelomereResult",
]
