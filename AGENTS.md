# simulation-project — hard-square-yang-baxter-obstruction

Leer `.diagnostika/structure.json` y el snapshot `skills/` pertinente.
No mover raíces sin ADR en docs/decisions ni modificar snapshots sin actualizar lock.
Registrar claims y evidencia reales; no inventar resultados ni aprobaciones.
La revisión semántica se activa en documentos argumentativos y manuscritos.
Lean y revisión semántica son garantías independientes.
Usar Python >=3.12 y uv; ejecutar validador y pytest antes de entregar.
No publicar fuentes privadas, datos ni credenciales. Sin migraciones automáticas.

## Reglas propias de este repositorio

- `verify_certificates.py` es la única fuente del código y es idéntico byte a byte
  al archivo auxiliar de arXiv (docs/decisions/0001-script-canonico-en-raiz.md).
- Al cambiar una comprobación, actualizar en el mismo commit `CLAIMS_MAP.md`,
  `claims/claim_registry.yaml` y el número esperado en `config/simulations/scenarios.yaml`.
- Una comprobación que pasa respalda un cálculo, no un teorema: las pruebas están
  en el paper.
- Las versiones siguen al paper (`v1.0.0` = arXiv v1). Cada release se archiva en
  Zenodo; no mover tags publicados.
