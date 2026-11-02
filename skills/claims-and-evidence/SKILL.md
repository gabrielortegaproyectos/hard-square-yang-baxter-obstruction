---
name: claims-and-evidence
description: Registrar afirmaciones científicas y asociar evidencia trazable al redactar resultados, tablas o conclusiones.
---

# claims-and-evidence

Una afirmación comunicable recibe un claim_id estable y único en
claims/claim_registry.yaml. Describe statement, claim_type, scope, status,
limitations, evidence y provenance. Tipos: finding, methodology, interpretation,
hypothesis, proposal, open_question. Estado científico: proposed, supported,
refuted, inconclusive. lifecycle independiente: active, historical, superseded.
No confundas una hipótesis con un resultado observado ni supported con verdad universal.

Asigna evidence_id a cada tabla, figura, manifest, run, config, documento,
presentación o resultado versionado en claims/evidence_map.yaml. Relaciona IDs
desde el claim. Cada evidencia registra path/URI, kind, run_id cuando aplique,
config, commit y checksum del artefacto cuando sea verificable. Usa ubicaciones
accesibles a los lectores autorizados; nunca publiques credenciales ni datos privados.
Un documento cita claim_id/evidence_id y localizador de tabla/página/slide.

Antes de promover supported, abre la evidencia, confirma que respalda el alcance,
anota incertidumbre y limitaciones, y comprueba que el número se reproduce desde
run/config. Una lista de rutas no constituye por sí sola una prueba científica.
No sobrescribas evidencia histórica: agrega nuevo ID/run y supersedes al reemplazar.
Lee [references/example.yaml](references/example.yaml) al registrar el primer claim;
es un ejemplo sintético, no un resultado listo para afirmar.

El helper `scripts/validate_references.py` comprueba IDs, referencias y contratos
de investigación aplicada cuando existen. No aprueba soporte científico ni prosa.
