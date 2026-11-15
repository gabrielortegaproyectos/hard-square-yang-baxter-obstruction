# simulation-project — hard-square-yang-baxter-obstruction

Leer `.diagnostika/structure.json` y el snapshot `skills/` pertinente.
No mover raíces sin ADR en docs/decisions ni modificar snapshots sin actualizar lock.
Registrar claims y evidencia reales; no inventar resultados ni aprobaciones.
La revisión semántica se activa en documentos argumentativos y manuscritos.
Lean y revisión semántica son garantías independientes.
Usar Python >=3.12 y uv; ejecutar validador y pytest antes de entregar.
No publicar fuentes privadas, datos ni credenciales. Sin migraciones automáticas.

## Reglas propias de este repositorio

- `verify_certificates.py` es la única fuente de las comprobaciones de arXiv v1 y es
  idéntico byte a byte a su archivo auxiliar (docs/decisions/0001-script-canonico-en-raiz.md).
  Los cálculos posteriores son módulos de `src/project/`, cada uno con su protocolo en
  `config/simulations/` y su informe en `docs/simulations/`.
- Al cambiar una comprobación, actualizar en el mismo commit `CLAIMS_MAP.md`,
  `claims/claim_registry.yaml` y lo esperado en su protocolo de `config/simulations/`
  (`scenarios.yaml` para `verify_certificates.py`).
- Una comprobación que pasa respalda un cálculo. Si el paper prueba el enunciado a mano,
  la prueba está en el paper; si lo presenta como calculado, el cálculo es su respaldo
  y su informe en `docs/simulations/` da el método y los límites.
- Las versiones siguen al paper (`v1.0.0` = arXiv v1). Cada release se archiva en
  Zenodo; no mover tags publicados.
