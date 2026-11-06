# Matrices de transferencia que conmutan con la de hard squares, para N pequeño

Issue #1. Protocolo: [`config/simulations/commuting-small-n.yaml`](../../config/simulations/commuting-small-n.yaml).
Código: [`src/project/commuting.py`](../../src/project/commuting.py) (`uv run python -m project.commuting`,
con `--long` para N = 7 y el tangente hasta N = 12). Pruebas: [`tests/test_commuting.py`](../../tests/test_commuting.py).

Es un cálculo exploratorio y no está en el archivo auxiliar de arXiv;
`verify_certificates.py` no cambia. El autor pidió llevarlo a la Sección 5.2 de la versión
para J. Phys. A (PR #54 del repositorio del paper), que sigue pendiente de su lectura.

## Pregunta

La Sección 5.2 de la versión para J. Phys. A (PR #54 del repositorio del paper) propone
resolver, para N pequeño,

    [T_N(Q), T_N(W)] = 0,

donde Q es una cara hard-square pura y W una cara binaria con pesos en las siete caras
admisibles. Es la prueba que sugiere Lochak–Maillard (J. Math. Phys. 27 (1986) 593): no
supone ninguna relación local. El Teorema B supone una ecuación de Yang–Baxter dentro
de una familia analítica; esta prueba no.

## Convenciones

Son las de `verify_certificates.py`. La cara (a b ; c d) tiene fila superior a b y fila
inferior c d, y los siete pesos son w0..w6 = 0000, 0001, 0010, 0100, 0110, 1000, 1001. La
cara pura es Q = (1, D, C, B, BC, A, AD), con z = ABCD ≠ 0. En filas hard-square
periódicas de largo N, T_N(W)(x, y) = ∏_i W[y_i y_{i+1} ; x_i x_{i+1}], con x la fila
inferior.

## Reducción

- **Solo importa z.** T_N(Q) = D_CD C_N D_AB, con C_N la matriz de compatibilidad de
  filas hard-square y D_t = diag(t^|x|). Conjugar por D_CD lleva el problema para Q al
  problema para Q0 = (1,1,1,1,1,z,z), con W reemplazada por W (CD)^(a−c).
- **Simetrías triviales.** La escala multiplica T_N por λ^N. Los dos gauges de esquina
  W ↦ W t^(a−b) s^(c−d) no cambian T_N. Son tres direcciones que siempre están en el
  conjunto de soluciones.
- **Soluciones evidentes.** δ_{b,c} = (1,1,0,0,1,1,1) da el desplazamiento S y
  δ_{a,d} = (1,0,1,1,1,0,1) da S⁻¹. Ambas conmutan con T_N(Q).

El módulo comprueba estas identidades en forma simbólica para N = 2, …, 6.

## Resultados

Para cada N, la tabla da el número de filas, el número de ecuaciones distintas (una por
par de filas salvo rotación) y las soluciones con T_N(W) ≠ 0, salvo escala y gauges, con
z simbólico y z ≠ 0. La última columna es la dimensión del espacio tangente en Q0, con las
tres direcciones triviales incluidas.

| N | filas | ecuaciones | soluciones | tangente en Q0 |
|---|---|---|---|---|
| 2 | 3 | 4 | componentes de dimensión 2 en soporte completo, y otras | 5 |
| 3 | 4 | 4 | componentes de dimensión 2 | 5 |
| 4 | 7 | 13 | solo Q0, δ_{b,c}, δ_{a,d} como conjunto; Q0 es un punto doble | 4 |
| 5 | 11 | 25 | exactamente Q0, δ_{b,c}, δ_{a,d}, puntos reducidos | 3 |
| 6 | 18 | 58 | lo mismo | 3 |
| 7 | 29 | 112 | lo mismo (`--long`) | 3 |

1. **Desde N = 4 no hay compañeros nuevos.** Para N = 4, 5, 6, 7 y todo z ≠ 0, complejo
   incluido, sea W una cara binaria de siete pesos con T_N(W) ≠ 0. Entonces T_N(W)
   conmuta con T_N(Q) si y solo si es múltiplo de T_N(Q), de S o de S⁻¹.
   - El cálculo recorre los 127 soportes no vacíos; los que dan T_N(W) = 0 quedan fuera.
     En cada uno es una base de Gröbner sobre Q[w, z] saturada por z y por los pesos
     libres, así que vale para todo z ≠ 0 a la vez: no aparece ningún valor especial de z,
     tampoco z = −1.
   - Con N = 4, las caras (1,1,0,0,1,1,−1) y (1,0,1,1,−1,0,1) también son soluciones, pero
     dan las mismas matrices T_4 = S y T_4 = S⁻¹. Desde N = 5 ya no son soluciones.
2. **Q es rígida al primer orden desde N = 5.** Con N = 5 y N = 6, los factores
   invariantes sobre Q[z] de las ecuaciones tangentes en las direcciones no triviales son
   1, 1, z², z². Así, el tangente es trivial para todo z ≠ 0. Con N = 4 hay una dirección
   tangente de más, y la base de Gröbner muestra que no se integra: el punto es doble.
   Con aritmética racional exacta, el tangente es trivial para 5 ≤ N ≤ 12 en
   z = 2, −1, 1/3, −1/4 y 7/5.
3. **Para N = 2 y N = 3 las condiciones son débiles.** Las soluciones forman
   componentes de dimensión 2, salvo simetrías. Son las condiciones algebraicas del tipo
   de Lochak–Maillard para tamaños chicos, y desde N = 4 dejan solo las soluciones
   evidentes.

### Control positivo: hard hexagons

Con la cara de hard hexagons (1,1,1,1,1,z,0) en lugar de Q0, a z = 2, el mismo código
encuentra en el soporte completo una curva de compañeros que conmutan.

- La curva es la misma para N = 4, 5 y 6.
- Su ideal reduce a cero el numerador de F(z_W, X, Y) = z_W XY(XY − X − Y) − (X − 1)(Y − 1),
  con X = w0w6/(w1w5), Y = w0w4/(w2w3) y z_W = w1w2w3w5/w0⁴. Es decir, está sobre la
  superficie integrable de Baxter.
- La actividad z_W varía a lo largo de la curva, como en las familias que conmutan de
  Baxter.

El método, entonces, encuentra una familia que conmuta cuando existe, y su ausencia para
hard squares no es un artefacto del código. Además, las dos partes del cálculo concuerdan:
el tangente (código lineal) y la variedad completa (Gröbner) dan un punto doble con
N = 4 y un punto reducido desde N = 5.

### Cruce independiente

[`src/project/commuting_crosscheck.py`](../../src/project/commuting_crosscheck.py)
(`uv run python -m project.commuting_crosscheck`, unos 40 s) repite el cálculo de siete
pesos sin importar nada del primer módulo:

- **Otras ecuaciones.** Usa todas las entradas del producto explícito de matrices
  T_N(Q0) T_N(W) − T_N(W) T_N(Q0), sin reducir por rotación, con una enumeración de filas
  escrita de nuevo.
- **Otra resolución.** En cada soporte fija a 1 otro conjunto de pesos. Usa bases grevlex
  con el algoritmo F5B, sin eliminación. En cada punto clasifica T_N(W) comparándola con
  T_N(Q0), S y S⁻¹, y lee la multiplicidad de la dimensión del cociente.
- **Otro tangente.** Calcula el jacobiano de las ecuaciones explícitas por derivación
  simbólica.
- **Otras actividades.** Usa z = 3, −2/7 y 5/2 en lugar de z simbólico; con z simbólico,
  F5B es demasiado lento.

Para N = 4, 5 y 6 da lo mismo que el primer método:
- solo las tres clases evidentes;
- Q0 doble con N = 4 y reducido desde N = 5;
- el mismo rango del jacobiano.

Con N = 4, cada soporte de desplazamiento tiene dos órbitas, δ y su variante con signo.
Ambas dan T_4 = S (o S⁻¹).

## El alfabeto A_4, cerca de hard squares

[`src/project/commuting_a4.py`](../../src/project/commuting_a4.py)
(`uv run python -m project.commuting_a4`, alrededor de un minuto; `--long` agrega M = 10).

En A_4 una cara es un par (W^0, W^1) de caras binarias, una por tipo: catorce pesos. Las
filas viven en cilindros de circunferencia par M, y cada fila lleva la paridad de su primera
altura; filas consecutivas tienen paridades opuestas.

- **Reducción.** Factores de esquina y una semejanza diagonal de las filas, ambos realizados
  en las caras, llevan toda presentación de actividad z a Q0 ⊕ Q0, la misma cara Q0 en los
  dos tipos. Está comprobado en forma simbólica.
- **Simetrías triviales.** Son seis:
  - la escala;
  - W^0 ↦ μW^0 con W^1 ↦ W^1/μ;
  - cuatro gauges de esquina, dos en la fila de arriba y dos en la de abajo.

  Todas multiplican T_M por un escalar.
- **Pregunta local.** Si el tangente en Q del conjunto de soluciones tiene dimensión 6, cerca
  de Q las soluciones son la órbita de las simetrías. Entonces toda familia analítica que
  pase por Q y conmute tiene T_M proporcional a T_M(Q). Es la versión con matrices de
  transferencia del Teorema B en A_4.
- **Por qué solo local.** Con catorce pesos, la eliminación global del caso binario no es
  factible con sympy. No se clasifican los compañeros lejos de Q.

| M | filas | tangente en Q |
|---|---|---|
| 4 | 14 | 7: una dirección más que las triviales |
| 6 | 36 | 6 para todo z ≠ 0 (el mcd de los menores 8 × 8 divide a z¹²) |
| 8 | 94 | 6 en z = 2, −1, 1/3, −1/4, 7/5 |
| 10 | 246 | 6 en los mismos cinco valores (`--long`) |

- **Sin usar la reducción.** Una presentación escalonada, con factores de esquina distintos
  en cada tipo y s₁/s₀ = 3/2, da lo mismo con M = 6.
- **Control.** Con hard hexagons en ambos tipos, el tangente tiene dimensión 7 (M = 6 y 8).
  La dirección extra es no escalonada: es la de la familia de Baxter en el alfabeto binario.
- **M = 4 es chico.** Pasa lo mismo que con N = 4 en el caso binario.

## Qué se sigue

- **No hay familias que pasen por Q.** Sea W(u) una familia analítica de caras binarias de
  siete pesos con W(u₀) = Q. Si T_N(W(u)) conmuta con T_N(Q) para un solo N entre 4 y 7,
  entonces T_N(W(u)) es múltiplo de T_N(Q) para todo u. Cerca de Q, todos los pesos son
  no nulos y el conjunto de soluciones es la órbita de Q. No se supone Yang–Baxter, ni
  conmutación mutua dentro de la familia, ni N arbitrario.
- **Tampoco hay compañeros aislados.** Ninguna cara de siete pesos, cerca de Q o lejos,
  conmuta con hard squares para esos tamaños, salvo las evidentes. Esto responde la
  pregunta de la Sección 5.2 en el caso de "ausencia".
- **En A_4, lo mismo cerca de hard squares.** Ninguna familia analítica de caras de A_4 que
  pase por una presentación tiene matrices de transferencia que conmuten con la de hard
  squares para M = 6 u 8, salvo las proporcionales.
- **Es complementario al Teorema B.** El teorema cubre familias de Yang–Baxter en los
  alfabetos binario y A_4. Este cálculo no supone relación local: en el binario cubre todas
  las caras de siete pesos, y en A_4 las familias que pasan por hard squares. Juntos dan la
  versión sin Yang–Baxter del Teorema B en los dos alfabetos del Teorema A, para estos
  tamaños.

## Límites

- **Qué cubre.**
  - Binario: siete pesos, N ≤ 7 y todo z ≠ 0.
  - A_4: solo cerca de una presentación, con M ≤ 8 (10 con `--long`).
  - No cubre alturas ocultas ni bloques.
  - Tampoco cubre el binario con 1∼1: ahí T_N(Q̂) = T_N(Q) ⊕ 0, y la conmutación no
    restringe las filas con vecinos ocupados. Ese caso, junto con el de nueve pesos, se
    calculó y se retiró porque no responde una pregunta del paper (commit b897c79).
- **Tamaño finito.** En el binario, el tangente vale para todo z solo con N = 5 y 6; para
  N ≤ 12 vale en cinco valores de z. En A_4 vale para todo z con M = 6.
- **Herramienta.** Todas las bases de Gröbner son de sympy. El cruce independiente usa otro
  algoritmo (F5B), otro orden y otro corte de gauge, pero el mismo sistema de álgebra
  computacional.
- **No es un teorema del paper.** Nada aquí afirma no-integrabilidad absoluta, ni dice
  nada sobre κ, la D-finitud o la frontera natural.

## Posibles extensiones

1. A_4 lejos de hard squares (eliminación global), con otro sistema de álgebra
   computacional.
2. Conmutación mutua dentro de una familia candidata, en lugar de conmutación con Q.
3. El grupo de inversión de Bellon–Maillard–Viallet (Sección 5.3).
