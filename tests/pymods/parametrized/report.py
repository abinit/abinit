"""
Comparison report of a parametrized template: comparison_report.json, the machine-readable result of
`run_params.py run/compare` (plan section 6). The HTML pages are written by htmlreport.py,
the same pages as runtests.py.

Standard library only.
"""
from __future__ import annotations

import datetime
import json
import os
from typing import Any

from .execute import STATUSES, InstanceOutcome


def build_report(template: str, reference: str | None, workdir_root: str, command: str,
                 outcomes: list[InstanceOutcome]) -> dict[str, Any]:
    counts = {s: sum(1 for o in outcomes if o.status == s) for s in STATUSES}
    return {
        "template": os.path.abspath(template),
        "reference": os.path.abspath(reference) if reference else None,
        "workdir_root": os.path.abspath(workdir_root),
        "command": command,
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "summary": counts,
        "runs": [o.as_dict() for o in outcomes],
    }


def write_json(report: dict[str, Any], path: str) -> None:
    with open(path, "w") as fh:
        json.dump(report, fh, indent=2)


def write_reports(report: dict[str, Any], json_path: str) -> str:
    """Write the JSON report. Return its path."""
    os.makedirs(os.path.dirname(os.path.abspath(json_path)), exist_ok=True)
    write_json(report, json_path)
    return json_path
