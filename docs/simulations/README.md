# Simulaciones
Cada escenario fija modelo/version, parámetros, distribuciones, semillas, configuración, entorno, commit y salida. Registrar run_id y checksums. Sensibilidad compara escenarios con criterios comunes; no asumir validez por reproducibilidad.

## En este proyecto

Los escenarios de `config/simulations/scenarios.yaml` son corridas deterministas de
verificación exacta, no simulaciones estocásticas: no hay semillas ni distribuciones.
El entorno lo fija `uv.lock`, y la salida de cada corrida de CI se publica como
artefacto `verify-output` (contrato en `contracts/simulation-runs.yaml`). El análisis de
sensibilidad no aplica: una comprobación exacta pasa o falla.
