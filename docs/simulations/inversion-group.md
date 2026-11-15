# El grupo de inversión de las caras binarias

Issue #3. Protocolo: [`config/simulations/inversion-group.yaml`](../../config/simulations/inversion-group.yaml).
Código: [`src/project/inversion_group.py`](../../src/project/inversion_group.py) (`uv run python -m project.inversion_group`).
Pruebas: [`tests/test_inversion_group.py`](../../tests/test_inversion_group.py).

Acompaña a la Proposición 5.4 y a la Observación 5.5 de la versión para J. Phys. A del
paper (PR #54 del repositorio del paper), que las prueba a mano. No está en el archivo
auxiliar de arXiv; `verify_certificates.py` no cambia.

## Pregunta

Para modelos de vértices y de espines, Bellon, Maillard y Viallet mostraron que las dos
relaciones de inversión generan un grupo infinito de transformaciones birracionales, cuyas
órbitas localizan las variedades donde vive el parámetro espectral (Phys. Lett. B 260
(1991) 87; Phys. Lett. A 157 (1991) 343). En un modelo de caras, una relación de inversión
invierte el peso leído como matriz en las dos alturas de una diagonal, con las de la otra
fijas (Maillard y Garel, J. Phys. A 17 (1984) 1257, Sección 2), y hay una por diagonal.

La pregunta es qué grupo generan sobre las caras binarias no nulas en las siete caras
admisibles, es decir, sobre el modelo de hard squares con interacciones diagonales, y qué
pasa en A_4.

## Convenciones

Son las de `verify_certificates.py`. La cara (a b ; c d) tiene fila superior a b y fila
inferior c d, y los siete pesos son w0..w6 = 0000, 0001, 0010, 0100, 0110, 1000, 1001. Las
coordenadas son

    z = w1 w2 w3 w5 / w0^4,    X = w0 w6 / (w1 w5),    Y = w0 w4 / (w2 w3),
    xi = 1/X,    eta = 1/Y,

que determinan W salvo un escalar y factores de esquina. La superficie de Baxter es
F(z, X, Y) = z X Y (X Y - X - Y) - (X - 1)(Y - 1) = 0. La cara hard-square pura es
Q = (1, D, C, B, BC, A, AD), con X = Y = 1.

## Qué se comprueba

Todas las comprobaciones son igualdades de funciones racionales en pesos simbólicos.

1. **Las matrices.** Con (a, d) fijas, la cara leída en (b, c) da la matriz
   [[w0, w2], [w3, w4]] y los números w1, w5 y w6. Con (b, c) fijas, leída en (a, d), da
   [[w0, w1], [w5, w6]] y los números w2, w3 y w4. El módulo las construye desde la lista
   de caras admisibles.
2. **Las inversiones.** I_1 invierte las primeras e I_2 las segundas. Con
   om = w0 w4 - w2 w3, I_1(W) = (w4/om, 1/w1, -w2/om, -w3/om, w0/om, 1/w5, 1/w6), e I_2
   es análoga. Los determinantes son w2 w3 (Y - 1) y w1 w5 (X - 1): I_1 no existe en Y = 1,
   I_2 no existe en X = 1, y en una cara hard-square pura no existe ninguna de las dos.
3. **Acción en coordenadas.**

       I_1: (z, xi, eta) -> (eta^2 (1 - eta)^2 / z, (1 - eta)/xi, eta)
       I_2: (z, xi, eta) -> (xi^2 (1 - xi)^2 / z, xi, (1 - xi)/eta)

   Ambas son involuciones, sobre los pesos y sobre las coordenadas.
4. **Período cinco.** La sucesión zeta_0 = xi, zeta_1 = eta,
   zeta_{n+1} = (1 - zeta_n)/zeta_{n-1} cumple zeta_5 = xi y zeta_6 = eta. La composición
   I_2 I_1 lleva (zeta_0, zeta_1) a (zeta_2, zeta_3) y multiplica z por
   gamma(zeta_2)/gamma(zeta_1), con gamma(s) = s^2 (1 - s)^2.
5. **El grupo.** (I_2 I_1)^5 es la identidad sobre (z, xi, eta) y ninguna potencia menor
   lo es. La clausura bajo I_1 e I_2 tiene exactamente 10 aplicaciones: un grupo diédrico
   de orden 10. Sobre los siete pesos, (I_2 I_1)^5 W es W por un escalar y factores de
   esquina de producto 1.
6. **Superficie de Baxter.** F o I_1 = F/(z X^2 Y) y F o I_2 = F/(z X Y^2), así que las dos
   inversiones conservan F = 0.
7. **A_4.** Con un peso binario por tipo, las inversiones actúan tipo por tipo sobre xi_t y
   eta_t. La actividad de una partícula en la subred de paridad p es
   a_p = w5^(p) w1^(p) w3^(1-p) w2^(1-p) / (w0^(p) w0^(1-p))^2, y

       I_1: a_p -> (1 - eta_p)^2 eta_{1-p}^2 / a_p
       I_2: a_p -> (1 - xi_{1-p})^2 xi_p^2 / a_p

   Sea Pi_t = (1 - xi_t)(1 - eta_t)(xi_t + eta_t - 1)/(xi_t eta_t), el producto de los cinco
   términos de un período; el producto de los 1 - zeta_n es Pi_t^2. Entonces (I_2 I_1)^5 fija
   xi_t y eta_t y lleva (a_0, a_1) a (mu a_0, a_1/mu), con mu = (Pi_1/Pi_0)^2. Si los dos
   tipos tienen los mismos xi y eta, mu = 1.

La corrida imprime catorce líneas `[PASS]` y tarda unos 10 segundos.

## Qué respalda y qué no

- **Respalda** las fórmulas que usa la prueba de la Proposición 5.4 y el cálculo indicado
  en la Observación 5.5. La consecuencia que saca el paper es que, en este patrón de pesos,
  las órbitas del grupo de inversión son finitas: no producen las familias de Baxter ni
  descartan ninguna cara, y lo que el grupo distingue es la recta de hard squares puro,
  donde ninguna inversión existe.
- **No es la prueba.** Una comprobación que pasa respalda un cálculo, no un teorema; la
  prueba está en el paper.
- **No dice nada de otros patrones**, como las dieciséis caras binarias o la forma de
  vértices, donde el grupo puede ser infinito.
- **No afirma novedad.** No se ha hecho una revisión de literatura sobre este grupo.
- **Un solo sistema.** Todo se calcula con sympy, sin cruce con otro sistema de álgebra
  computacional. La prueba `test_exact_rational_point` sigue un punto racional fijo por los
  diez elementos del grupo, como comprobación aritmética aparte.
