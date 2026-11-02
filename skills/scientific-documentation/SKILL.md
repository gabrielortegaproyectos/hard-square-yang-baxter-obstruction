---
name: scientific-documentation
description: Escribir o revisar documentación científica, manuscritos e informes con resultados trazables y límites explícitos.
---

# scientific-documentation

Expón pregunta y alcance, metodología, resultados observados, interpretación,
limitaciones, hipótesis y preguntas abiertas como categorías distinguibles.
Cada cifra debe citar tabla/figura/artefacto con run/config y revisión de código;
cada conclusión científica debe referir un claim estable. No infieras causalidad
de una asociación ni generalización de una partición sin justificarla.

Escribe narrativa humana en docs/reporting/. Las tablas machine-readable y manifests
van en la capa de datos de resultados, no en docs ni en otra raíz reporting/.
Conserva versiones previas y marca active/historical/superseded con enlace al sucesor.
Al revisar un informe, comprueba unidades, denominadores, incertidumbre, población,
fecha de corte y consistencia entre texto, tablas y claims. Registra faltantes sin
inventar resultados. Para manuscritos consulta [references/manuscripts.md](references/manuscripts.md) y activa la skill `research` al organizar claims, párrafos, secciones e inferencias del manuscrito (también cuando el manuscrito viva en un proyecto ML o de matemática formal).
