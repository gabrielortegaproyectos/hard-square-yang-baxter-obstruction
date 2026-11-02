---
name: reproducibility
description: Preparar o verificar ejecuciones repetibles, provenance y entornos al producir resultados o entregas.
---

# reproducibility

Usa Python >=3.12 y uv; versiona pyproject.toml y uv.lock, nunca .venv ni .env.
Documenta comandos, versión de runtime, seed y límites de determinismo por hardware.
Para cada run inmutable conserva run_id, timestamps UTC, git SHA (y dirty si aplica),
checksum de entradas/contratos, config congelada, entorno/lock y rutas de salidas.
En ML guarda manifest en data/08_reporting/<run_id>/; ETL en data/09_quality_reports/<run_id>/.
Los manifests son datos de provenance, no una evaluación Criterium.

Prueba primero con datos sintéticos. Para entradas privadas explica recepción y
verificación de hashes sin publicarlas. Registra metadata/snapshots sólo cuando
provenance atraviese productos/repositorios y la extensión esté justificada.
Para dependencias incompatibles usa envs/<componente>/ con lock propio y un comando
documentado; evita hacerlo obligatorio. Reejecuta en otro destino/run_id, compara
métricas con tolerancia explícita y nunca sobrescribas el run original.
