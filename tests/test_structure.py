from pathlib import Path
import runpy
import shutil

ROOT = Path(__file__).resolve().parents[1]
validate = runpy.run_path(str(ROOT / 'scripts/validate_structure.py'))['validate']


def test_template():
    assert validate(ROOT) == []


def test_corruption_and_missing_contract(tmp_path):
    target = tmp_path / 'project'
    shutil.copytree(ROOT, target, ignore=shutil.ignore_patterns('.git', '.venv', '__pycache__', '.pytest_cache'))
    skill = next((target / 'skills').rglob('SKILL.md'))
    skill.write_text('drift', encoding='utf-8')
    assert any('snapshot' in e for e in validate(target))
    (target / '.diagnostika/structure.json').unlink()
    assert validate(target)
