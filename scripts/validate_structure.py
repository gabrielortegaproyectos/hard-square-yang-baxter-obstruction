"""Standalone structure/hash/reference checks, independent of the CLI."""
import hashlib
import json
from pathlib import Path, PurePosixPath
import runpy
import yaml


def validate(root):
    errors = []
    try:
        spec = json.loads((root / '.diagnostika/structure.json').read_text(encoding='utf-8'))
        meta = yaml.safe_load((root / '.diagnostika/project.yaml').read_text(encoding='utf-8'))
        lock = json.loads((root / '.diagnostika/skills.lock.json').read_text(encoding='utf-8'))
        if spec['project_type'] != meta['project']['type']:
            errors.append('Project identity mismatch')
        for raw in spec['required_paths']:
            path = PurePosixPath(raw)
            if path.is_absolute() or '..' in path.parts or ':' in raw or chr(92) in raw:
                errors.append('Unsafe path: ' + raw)
            elif not (root / raw).exists():
                errors.append('Missing: ' + raw)
        for name in spec['forbidden_roots']:
            if (root / name).exists():
                errors.append('Forbidden root: ' + name)
        if set(spec['skills']) != set(lock['skills']):
            errors.append('Skills mismatch')
        actual = {p.relative_to(root / 'skills').as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in (root / 'skills').rglob('*') if p.is_file()}
        if actual != lock['files']:
            errors.append('Skills snapshot differs from lock')
        if any(p.is_symlink() for p in root.rglob('*') if p.relative_to(root).parts[0] not in {'.venv', '.git', '.pytest_cache'}):
            errors.append('Links are unsupported')
        if (root / 'claims/claim_registry.yaml').exists():
            checker = root / 'skills/claims-and-evidence/scripts/validate_references.py'
            errors.extend(runpy.run_path(str(checker))['validate'](root))
        # Formal projects preserve their own registry and independent checks.
        formal = root / 'scripts/check_formalization.py'
        if formal.exists():
            import subprocess, sys
            result = subprocess.run([sys.executable, str(formal)], capture_output=True, text=True)
            if result.returncode:
                errors.append('Formalization: ' + result.stdout + result.stderr)
    except (OSError, ValueError, KeyError, TypeError, yaml.YAMLError) as error:
        errors.append('Invalid contract: ' + str(error))
    return errors


if __name__ == '__main__':
    errors = validate(Path(__file__).resolve().parents[1])
    print('\n'.join(errors) if errors else 'Structure and references valid; contextual review remains required')
    raise SystemExit(bool(errors))
