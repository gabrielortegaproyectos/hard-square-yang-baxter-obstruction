"""Run the verification script as a test: every check must pass."""
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "verify_certificates.py"


def run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True, check=False)


def test_default_checks_pass():
    result = run()
    assert result.returncode == 0, result.stdout + result.stderr
    assert "All checks passed." in result.stdout
    assert "[FAIL]" not in result.stdout


@pytest.mark.slow
def test_groebner_and_classification_pass():
    result = run("--groebner", "--classify")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "All checks passed." in result.stdout
