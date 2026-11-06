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

## Más allá de los siete pesos

[`src/project/commuting_extended.py`](../../src/project/commuting_extended.py)
(`uv run python -m project.commuting_extended`, unos 70 s; `--long` agrega z = −1/3).

### Nueve pesos: las caras de par vertical

Las caras (10;10) y (01;01) conservan las filas hard-square, pero permiten un sitio ocupado
en dos filas seguidas. Cada par vertical usa una de cada una, así que entran en T_N solo
por el producto w7·w8. Con ellas aparecen tres compañeros evidentes más:
- la identidad I, que viene de W = 1 en 0000, 1010 y 0101;
- K·S y K·S⁻¹: hard squares con la exclusión vertical cambiada por una diagonal, que son
  productos de matrices que conmutan.

Resultado a z = 2 (y a z = −1/3 con `--long`):

| N | clases de T_N(W) que conmutan con K = T_N(Q0) |
|---|---|
| 4 | K, S, S⁻¹, I, K·S, K·S⁻¹ y una solución esporádica de rango 1 |
| 5 | las mismas seis y dos esporádicas, de rango 8 y 11 |
| 6 | exactamente K, S, S⁻¹, I, K·S, K·S⁻¹ |

Con N = 6, T_6(W) es múltiplo de K solo en soportes sin par vertical (w7·w8 = 0). Así,
tampoco con nueve pesos hay una familia que pase por Q: cerca de Q0 las seis clases son
puntos distintos del espacio proyectivo de matrices, y una curva continua de matrices
tiene que quedarse en la clase de K.

### Dieciséis caras: el bloque fuera de hard squares no se ve

Con pesos en caras con dos sitios ocupados vecinos en horizontal, T_N(W) actúa sobre
todas las filas binarias.

- **Bloques.** La cara de hard squares extendida por cero, Q̂, cumple
  T_N(Q̂) = T_N(Q) ⊕ 0 (filas hard-square primero).
- **Reducción a nueve pesos.** det C_N = ±1 para 2 ≤ N ≤ 12, así que T_N(Q) es invertible
  si z ≠ 0. Entonces [T_N(Q̂), T_N(W)] = 0 equivale a tres condiciones:
  - el bloque hard-square de T_N(W), que solo depende de los nueve pesos de arriba,
    conmuta con T_N(Q);
  - los dos bloques fuera de la diagonal se anulan;
  - el bloque de las otras filas queda libre.
- **Una familia desacoplada.** Por eso hay una familia trivial que conmuta:
  W_u = Q̂ + u·[1111] tiene T_N(W_u) = T_N(Q) ⊕ u^N E, con E la proyección sobre la fila
  totalmente ocupada. Comprobado para N = 2, …, 6.
  - Sus matrices conmutan con T_N(Q̂) y entre sí para todo N.
  - Por la Proposición 4.3 del paper no viene de una familia de Yang–Baxter.

Es un ejemplo explícito de por qué el teorema de Lochak–Maillard pide positividad: con
pesos nulos, la conmutación para todo N no implica una relación local. En las dieciséis
caras, entonces, la conmutación sola no puede excluir familias. El argumento local de la
Proposición 4.3 es necesario, y en el sector hard-square la conmutación se reduce al caso
de nueve pesos.

## Qué se sigue

- **No hay familias que pasen por Q.** Sea W(u) una familia analítica de caras binarias de
  siete pesos con W(u₀) = Q. Si T_N(W(u)) conmuta con T_N(Q) para un solo N entre 4 y 7,
  entonces T_N(W(u)) es múltiplo de T_N(Q) para todo u. Cerca de Q, todos los pesos son
  no nulos y el conjunto de soluciones es la órbita de Q. No se supone Yang–Baxter, ni
  conmutación mutua dentro de la familia, ni N arbitrario.
- **Tampoco hay compañeros aislados.** Ninguna cara de siete pesos, cerca de Q o lejos,
  conmuta con hard squares para esos tamaños, salvo las evidentes. Esto responde la
  pregunta de la Sección 5.2 en el caso de "ausencia".
- **Con nueve pesos, lo mismo.** Con N = 6 (z = 2 y z = −1/3), los únicos compañeros son
  K, S^{±1}, I y K·S^{±1}, y ninguna familia pasa por Q.
- **Con dieciséis caras, la conmutación no alcanza.** Hay familias desacopladas que
  conmutan sin ser de Yang–Baxter, así que ahí el resultado es el de la Proposición 4.3,
  no este cálculo.
- **Es complementario al Teorema B.** El teorema cubre familias de Yang–Baxter en los
  alfabetos binario (con las 16 caras) y A_4. Este cálculo no supone relación local, y
  cubre los siete pesos admisibles y, con N = 6, los nueve que conservan las filas
  hard-square.

## Límites

- **Qué cubre.** Siete pesos para N ≤ 7 y todo z ≠ 0. Nueve pesos para N ≤ 6, a z = 2 y
  z = −1/3 (`--long`).
  - Con dieciséis caras, la conmutación con Q̂ no restringe el bloque de filas que no son
    hard-square.
  - No cubre alturas ocultas, bloques ni A_4.
- **N finito.** El tangente vale para todo z solo con N = 5 y 6; para N ≤ 12 vale en cinco
  valores de z.
- **Herramienta.** Todas las bases de Gröbner son de sympy. El cruce independiente usa otro
  algoritmo (F5B), otro orden y otro corte de gauge, pero el mismo sistema de álgebra
  computacional.
- **No es un teorema del paper.** Nada aquí afirma no-integrabilidad absoluta, ni dice
  nada sobre κ, la D-finitud o la frontera natural.

## Posibles extensiones

1. Nueve pesos con z simbólico y N = 7.
2. Conmutación mutua dentro de una familia candidata, en lugar de conmutación con Q.
3. El grupo de inversión de Bellon–Maillard–Viallet (Sección 5.3).
4. A_4: la conmutación con la matriz escalonada T_1 T_0 de la Observación 2.11.
