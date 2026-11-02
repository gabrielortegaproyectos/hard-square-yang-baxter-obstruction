# 0001 — Script canónico en la raíz

Fecha: 2026-10-25.

**Contexto.** El template `simulation-project` ubica el modelo en
`src/project/simulation.py`. Las comprobaciones de Paper 1 deben ser, además,
el archivo auxiliar de arXiv: un único archivo autocontenido que un árbitro pueda
ejecutar sin instalar este paquete.

**Alternativas.**
1. Mover el código a `src/project/` y generar el archivo auxiliar en cada release.
2. Mantener `verify_certificates.py` en la raíz como única fuente, y que
   `src/project/simulation.py` solo lo ejecute por escenario.

**Decisión.** Alternativa 2. El repositorio del paper trae el archivo desde un tag
fijo de este repositorio y comprueba su SHA-256. Así el código citado y el
archivo auxiliar son el mismo objeto.

**Revisión.** Revisar si el script se divide en módulos o si arXiv deja de ser el
canal del archivo auxiliar.
