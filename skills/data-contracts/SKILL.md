---
name: data-contracts
description: Definir o cambiar interfaces machine-readable entre productores y consumidores de datasets, targets, splits e inferencia.
---

# data-contracts

Crea contratos YAML/JSON en contracts/ cuando exista una interfaz concreta; no
rellenes contratos vacíos para satisfacer un listado. Narrativa en docs/, políticas
ejecutables en config/. ML suele necesitar dataset.yaml, target.yaml, splits.yaml
e inference.yaml; ETL source.yaml, input.yaml y output.yaml.

Especifica schema_version, propietario, productor, consumidor, unidad/fila, claves,
tipos, nulabilidad, dominios, unidades, privacidad, compatibilidad y validaciones.
Para target agrega disponibilidad temporal y transformación; splits define unidad
de separación y exclusiones; inference define entradas, salidas y casos inválidos.
Incluye un ejemplo sintético válido y uno inválido y prueba ambos en el límite.
Ante cambio incompatible versiona la interfaz y acuerda consumidores; no reinterpretes
silenciosamente un campo existente. Conserva checksum/revisión de la versión consumida.
