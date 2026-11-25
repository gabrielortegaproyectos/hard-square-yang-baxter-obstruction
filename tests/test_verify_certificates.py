"""Run the canonical verification script and enforce the declared scenarios."""
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "verify_certificates.py"


def run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True, check=False)


def assert_scenario(result: subprocess.CompletedProcess, expected: int) -> None:
    assert result.returncode == 0, result.stdout + result.stderr
    assert "All checks passed." in result.stdout
    assert "[FAIL]" not in result.stdout
    passes = [line for line in result.stdout.splitlines() if "[PASS]" in line]
    assert len(passes) == expected, result.stdout
    assert sum("(float)" in line for line in passes) == 2, result.stdout


def test_default_checks_pass():
    assert_scenario(run(), 40)


@pytest.mark.slow
def test_groebner_and_classification_pass():
    assert_scenario(run("--groebner", "--classify"), 42)
