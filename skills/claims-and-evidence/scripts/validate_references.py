"""Reference integrity only; never a scientific or editorial score."""
from pathlib import Path
from datetime import date
import re
import yaml


def validate(root: Path) -> list[str]:
    errors = []
    try:
        def read(path):
            return yaml.safe_load((root / path).read_text(encoding='utf-8')) or {}
        registry = read('claims/claim_registry.yaml')
        evidence = read('claims/evidence_map.yaml')
        claims, items = registry['claims'], evidence['evidence']
        def ids(rows, key):
            result = [r[key] for r in rows]
            if any(not isinstance(x, str) or not x.strip() for x in result):
                raise ValueError(f'Invalid {key}')
            if len(set(result)) != len(result):
                errors.append(f'Duplicate {key}')
            return set(result)
        claim_ids, evidence_ids = ids(claims, 'claim_id'), ids(items, 'evidence_id')
        for claim in claims:
            for ref in claim.get('evidence', []):
                if ref not in evidence_ids:
                    errors.append(f'Unknown evidence: {ref}')
            for ref in claim.get('dependencies', []):
                if ref not in claim_ids:
                    errors.append(f'Unknown claim: {ref}')
        for item in items:
            if 'path' in item:
                path = item['path']
                target = root / path
                if not target.resolve().is_relative_to(root.resolve()) or not target.is_file():
                    errors.append(f'Invalid evidence path: {path}')
        questions_path = root / 'docs/problem/research-questions.yaml'
        if questions_path.exists():
            questions = ids(read('docs/problem/research-questions.yaml')['questions'], 'question_id')
            sources = ids(read('contracts/sources.yaml')['sources'], 'source_id')
            protocols = sorted((root / 'docs/research').glob('DR*/protocol.yaml'))
            dr_ids = set()
            required = {'id', 'question_id', 'scope', 'inclusion_reason', 'search_protocol',
                        'consultation_dates', 'selection_criteria', 'source_ids', 'quality',
                        'applicability', 'limitations', 'claim_ids', 'supporting_evidence',
                        'contradicting_evidence', 'unresolved_questions', 'related_research',
                        'proposal_implications', 'status'}
            records = []
            for file in protocols:
                dr = read(file.relative_to(root))
                if required - dr.keys():
                    errors.append(f'Incomplete DR: {file.parent.name}')
                    continue
                if not re.fullmatch(r'DR[0-9]{2,}', dr['id']) or dr['id'] != file.parent.name or dr['id'] in dr_ids:
                    errors.append(f'Invalid DR identity: {dr["id"]}')
                dr_ids.add(dr['id'])
                records.append(dr)
                if dr['status'] not in {'planned', 'searching', 'synthesized', 'reviewed'}:
                    errors.append(f'Invalid DR status: {dr["id"]}')
                for consulted in dr['consultation_dates']:
                    try:
                        date.fromisoformat(str(consulted))
                    except ValueError:
                        errors.append(f'Invalid consultation date: {dr["id"]}')
                if dr['status'] in {'synthesized', 'reviewed'}:
                    if not all(dr[field] for field in ('scope', 'inclusion_reason', 'search_protocol', 'consultation_dates', 'source_ids', 'quality', 'applicability')):
                        errors.append(f'Incomplete completed DR: {dr["id"]}')
                    for name in ('findings.md', 'limitations.md'):
                        if not (file.parent / name).is_file():
                            errors.append(f'Missing DR narrative: {dr["id"]}/{name}')
                if dr['question_id'] not in questions:
                    errors.append(f'Unknown question: {dr["question_id"]}')
                for field, known in [('source_ids', sources), ('claim_ids', claim_ids),
                                     ('supporting_evidence', evidence_ids), ('contradicting_evidence', evidence_ids)]:
                    if not isinstance(dr[field], list) or any(x not in known for x in dr[field]):
                        errors.append(f'Invalid DR {field}: {dr["id"]}')
            for dr in records:
                if any(x not in dr_ids for x in dr['related_research']):
                    errors.append(f'Unknown related DR: {dr["id"]}')
            for row in read('docs/proposal/recommendations.yaml')['recommendations']:
                refs = row['claim_ids']
                if not refs or any(x not in claim_ids for x in refs):
                    errors.append('Recommendation needs known claim_ids')
                for field in ('scientific_evidence', 'technical_feasibility', 'commercial_suitability', 'institutional_decision'):
                    if not row.get(field):
                        errors.append(f'Recommendation missing {field}')
    except (OSError, ValueError, KeyError, TypeError, AttributeError, yaml.YAMLError) as error:
        errors.append(f'Invalid reference contract: {error}')
    return errors
