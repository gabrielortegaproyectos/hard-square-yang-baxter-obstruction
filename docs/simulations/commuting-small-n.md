# Matrices de transferencia que conmutan con la de hard squares, para N pequeño

Issue #1. Protocolo: [`config/simulations/commuting-small-n.yaml`](../../config/simulations/commuting-small-n.yaml).
Código: [`src/project/commuting.py`](../../src/project/commuting.py) (`uv run python -m project.commuting`,
con `--long` para N = 7 y el tangente hasta N = 12). Pruebas: [`tests/test_commuting.py`](../../tests/test_commuting.py).

Es un cálculo exploratorio y no está en el paper ni en el archivo auxiliar de arXiv;
`verify_certificates.py` no cambia. Entra al manuscrito solo después de que el autor lo revise.

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

## Qué se sigue

- **No hay familias que pasen por Q.** Sea W(u) una familia analítica de caras binarias de
  siete pesos con W(u₀) = Q. Si T_N(W(u)) conmuta con T_N(Q) para un solo N entre 4 y 7,
  entonces T_N(W(u)) es múltiplo de T_N(Q) para todo u. Cerca de Q, todos los pesos son
  no nulos y el conjunto de soluciones es la órbita de Q. No se supone Yang–Baxter, ni
  conmutación mutua dentro de la familia, ni N arbitrario.
- **Tampoco hay compañeros aislados.** Ninguna cara de siete pesos, cerca de Q o lejos,
  conmuta con hard squares para esos tamaños, salvo las evidentes. Esto responde la
  pregunta de la Sección 5.2 en el caso de "ausencia".
- **Es complementario al Teorema B.** El teorema cubre familias de Yang–Baxter en los
  alfabetos binario (con las 16 caras) y A_4. Este cálculo no supone relación local, pero
  solo cubre los siete pesos admisibles.

## Límites

- **Solo siete pesos.** Hay caras binarias que no cubre:
  - Las caras de par vertical (10;10) y (01;01) conservan las filas hard-square, pero
    cada par vertical aporta dos de ellas. Por eso el tangente no las ve, y hace falta
    un cálculo global con nueve pesos.
  - Las caras con un par horizontal hacen actuar T_N(W) sobre todas las filas binarias.
    Ahí la matriz de Q extendida por cero es T_N(Q) ⊕ 0, y donde T_N(Q) es invertible la
    conmutación además anula los bloques fuera de la diagonal.

  Tampoco cubre alturas ocultas, bloques ni A_4.
- **N finito.** El cálculo global llega a N ≤ 7. El tangente vale para todo z solo con
  N = 5 y 6; para N ≤ 12 vale en cinco valores de z.
- **Herramienta.** Las bases de Gröbner son de sympy (lex). Las respaldan el control
  positivo y la concordancia entre el cálculo tangente y el global; no se cruzaron con
  otro sistema de álgebra computacional.
- **No es un teorema del paper.** Nada aquí afirma no-integrabilidad absoluta, ni dice
  nada sobre κ, la D-finitud o la frontera natural.

## Posibles extensiones

1. Nueve pesos (siete más las dos caras de par vertical) sobre filas hard-square, y las
   dieciséis caras con la condición de bloques.
2. Conmutación mutua dentro de una familia candidata, en lugar de conmutación con Q.
3. El grupo de inversión de Bellon–Maillard–Viallet (Sección 5.3).
