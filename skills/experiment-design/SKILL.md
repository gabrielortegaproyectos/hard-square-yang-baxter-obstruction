---
name: experiment-design
description: Diseñar o ejecutar experimentos ML o simulaciones con protocolo y decisiones congelados antes de mirar resultados.
---

# experiment-design

Antes de ejecutar fija pregunta, hipótesis, unidad estadística, target o resultado,
población, datos disponibles al predecir, baseline, split y métricas con dirección.
Registra presupuesto, seeds, criterios de selección/parada y análisis de sensibilidad.
Usa config/<protocol_id>.yaml con versión; los cambios posteriores crean otro ID
con motivo. No adaptes el protocolo al holdout después de observarlo.

En ML separa preparación, ajuste y evaluación; ajusta transformaciones sólo con train.
Guarda hiperparámetros y checkpoints bajo config y data/06_models/<run_id>.
Asocia outputs a un manifest y registra intentos fallidos para evitar selección sesgada.
Una corrida nueva no justifica otra raíz. Al terminar reporta observado, desviaciones,
incertidumbre y alcance; registra claims sólo con evidencia localizada.
