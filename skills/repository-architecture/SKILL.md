---
name: repository-architecture
description: Ubicar artefactos y decidir extensiones de estructura antes de agregar carpetas o experimentos.
---

# repository-architecture

Lee .diagnostika/project.yaml, .diagnostika/structure.json y docs/architecture.md.
Distribuye cada experimento bajo config/ (decisiones), src/ (lógica), scripts/
(entradas), data/<capa>/<run_id>/ (artefactos), docs/reporting/ (interpretación)
y claims/ (afirmaciones). No crees experiments-v2, pairwise ni new-analysis raíz.

Antes de agregar una raíz, documenta en docs/decisions/ el problema generalizable,
por qué no cabe en una raíz existente, semántica, dueño y política de privacidad.
Actualiza el contrato estructural y la documentación mediante un cambio revisable.
metadata/ es extensión para snapshots transversales de recepción: IDs, hashes,
versiones de contrato y commit upstream; no vuelques datos personales ni rutas secretas.
envs/ sólo se justifica con dependencias incompatibles: un pyproject y uv.lock por entorno.
No uses training/ para configs ni reporting/ para duplicar data/docs.
