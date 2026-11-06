# A local Yang–Baxter obstruction for the hard-square model — computations

This repository reproduces every finite computation in the paper
*A local Yang–Baxter obstruction for the hard-square model*
(G. Ortega and E. Peña-Ancavil,
[arXiv:2611.06780](https://arxiv.org/abs/2611.06780)).

Archived on Zenodo: version 1.0.0,
[doi:10.5281/zenodo.23892182](https://doi.org/10.5281/zenodo.23892182); all
versions, [doi:10.5281/zenodo.23892181](https://doi.org/10.5281/zenodo.23892181).

The paper shows that no analytic Yang–Baxter family of face weights passes
through hard squares at nonzero activity on either of the two connected height
alphabets that encode the model without hidden heights. The computations here
check:

- the classification of encodings;
- the polynomial certificates on the binary alphabet and on `A_4`;
- the fibre-sum presentations;
- the transfer-matrix identity;
- the blocked alphabet and the binary alphabet with `1~1`.

[`CLAIMS_MAP.md`](CLAIMS_MAP.md) lists which statement each check covers.

## Run

With [uv](https://docs.astral.sh/uv/), whose lockfile pins the environment:

```bash
uv run python verify_certificates.py
```

```bash
uv run python verify_certificates.py --groebner --classify
```

```bash
uv run pytest
```

With plain Python (the script itself runs on 3.10 or newer; the project pins 3.12), install `sympy` and `networkx` and
run `python verify_certificates.py`.

- **Default run:** about 15 seconds.
- **`--groebner --classify`:** about one minute. These two options add the
  Gröbner cross-check of the near-full-support theorem and the brute-force
  illustration of the classification.
- **Expected output:** 40 checks marked `[PASS]` in the default run, 42 with `--groebner --classify`, and `All checks passed.`
- **Arithmetic:** all checks are exact except the two calibration checks
  against the critical RSOS model, which are floating point and labelled as
  such.

## Exploratory computation: commuting transfer matrices

`src/project/commuting.py` asks which binary seven-weight faces have row transfer
matrices commuting with that of hard squares, `[T_N(Q), T_N(W)] = 0`, without assuming
any local relation (the question of Section 5.2 of the J. Phys. A version). For
`N = 4, …, 7` and every activity `z != 0`, only the evident ones do: multiples of `Q`
and the shift and its inverse. A positive control with hard hexagons recovers a curve of
commuting partners on Baxter's integrable surface. This computation is not part of the
paper or of the arXiv ancillary file. Method, results and limits are in
[`docs/simulations/commuting-small-n.md`](docs/simulations/commuting-small-n.md).

```bash
uv run python -m project.commuting
```

- **Default run:** about 20 seconds (`N <= 6`, tangent space up to `N = 9`).
- **`--long`:** about 3 minutes; adds `N = 7` and the tangent space up to `N = 12`.
- **Independent cross-check:** `uv run python -m project.commuting_crosscheck` (about 40 seconds).
- **The alphabet A_4:** `uv run python -m project.commuting_a4` (about a minute; `--long`
  adds circumference 10). Near a presentation of hard squares on A_4, only the trivial
  directions commute, so no analytic family of A_4 face weights through hard squares has
  commuting transfer matrices; hard hexagons has one more direction, Baxter's.

## Relation to the paper

- **Same script as on arXiv.** `verify_certificates.py` is the ancillary file
  of the arXiv submission, byte for byte (LF line endings).
- **Releases follow the paper versions:** `v1.0.0` for arXiv v1, `v1.1.0`
  for a revision, and so on. Each release is archived on Zenodo with its own
  DOI.

## Repository structure

The repository follows the Diagnostika Labs template `simulation-project`
(`v0.1.0-rc.2`, recorded in `.diagnostika/project.yaml`):

| Path | Role |
|---|---|
| `verify_certificates.py` | the canonical script (see `docs/decisions/0001-script-canonico-en-raiz.md`) |
| `config/simulations/scenarios.yaml` | the two scenarios: `default` (40 checks) and `full` (42 checks) |
| `contracts/simulation-runs.yaml` | what a recorded run must declare |
| `claims/` | the computational claims and the evidence that supports them |
| `src/project/simulation.py` | runs a scenario: `uv run python -m project.simulation full` |
| `src/project/commuting.py` | the exploratory commutation computation (protocol in `config/simulations/commuting-small-n.yaml`) |
| `skills/` | pinned agent-skill snapshot (hashes in `.diagnostika/skills.lock.json`) |

Validate the structure with `uv run python scripts/validate_structure.py`.

## Citation

See [`CITATION.cff`](CITATION.cff); GitHub shows it as "Cite this repository".

## License

MIT. See [`LICENSE`](LICENSE).
