"""Entry point for the verification scenarios of Paper 1.

The canonical code is the single-file script ``verify_certificates.py`` at the
repository root. It is byte-identical to the arXiv ancillary file and stays
self-contained for that reason (docs/decisions/0001-script-canonico-en-raiz.md).
This module only runs it with the arguments of a scenario declared in
``config/simulations/scenarios.yaml``.
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "verify_certificates.py"
SCENARIOS = ROOT / "config" / "simulations" / "scenarios.yaml"


def scenario_args(scenario_id: str) -> list[str]:
    """Return the command-line arguments of a declared scenario."""
    config = yaml.safe_load(SCENARIOS.read_text(encoding="utf-8"))
    for scenario in config["scenarios"]:
        if scenario["scenario_id"] == scenario_id:
            return list(scenario["arguments"])
    raise KeyError(f"Unknown scenario: {scenario_id}")


def run(scenario_id: str = "default") -> int:
    """Run one scenario; return the script's exit code (0 when all checks pass)."""
    saved = sys.argv
    sys.argv = [str(SCRIPT), *scenario_args(scenario_id)]
    try:
        runpy.run_path(str(SCRIPT), run_name="__main__")
    except SystemExit as stop:
        return int(stop.code or 0)
    finally:
        sys.argv = saved
    return 0


if __name__ == "__main__":
    raise SystemExit(run(sys.argv[1] if len(sys.argv) > 1 else "default"))
