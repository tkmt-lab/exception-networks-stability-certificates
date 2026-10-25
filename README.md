# Stability certificates for two exceptional reaction networks

Version 1.0.0.

This package contains the two reaction networks in Tang–Zhang's Example 4,
two derivations of their positive-equilibrium bounds, and exact certificates
for three fixed choices of rate constants. The certificates give rational
root-isolation intervals and Routh–Hurwitz signs. Two implementations recompute
them from the reaction table. These files can serve as test inputs for
steady-state and stability calculations.

For either network, with mass-action kinetics and **all six rate constants
strictly positive**:

- There are at most three distinct positive equilibria.
- If there are three distinct positive equilibria, exactly one is exponentially stable.
- There is at most one exponentially stable positive equilibrium, including when
  other positive equilibria are degenerate.

Here positive means that all three concentrations are strictly positive.
Exponential stability means that every eigenvalue of the Jacobian has negative
real part. These statements do not assert that three equilibria exist for every
choice of rates. The proofs are in [the ratio/y derivation](docs/ratio-y-derivation.md)
and [the z derivation](docs/z-derivation.md); [claims.json](claims.json) records
their assumptions and the separate finite-computation claims.

## Run the checks

With Python 3.14 and the versions in [requirements.txt](requirements.txt) installed,
run this command from the package directory:

```sh
python3 -B verify.py
```

For a new environment, first run:

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Dependency installation uses the network. The verification command itself runs
offline, with one checker subprocess at a time and a 120-second limit for each of its six
subprocesses. It requires a writable copy, creates temporary files beside the
runner, and removes them on completion. It leaves the saved certificates intact.
Exit status 0 means every comparison passed; a failed check returns nonzero.
Run without Python's `-O` option, because the archived programs use assertions.

The final result is:

```text
PASS: 2 networks, 3 samples, 9 equilibria; both implementations agree
```

The runner checks file hashes, constructs the reactant and stoichiometric matrices,
checks the symbolic identities, recomputes every supplied certificate, and matches
the equilibria between the two parametrizations. Timing and interpreter-version
fields in historical records are not required to match. Each implementation's
recomputed mathematical fields, including rational intervals, must match its own
saved certificates exactly. The reconciliation checks that the two parametrizations
identify the same equilibria and stability signs. The written arguments about
positive branches and root signs remain prose proofs; this command is not a proof
assistant checking every step. [VERIFICATION-RECORD.md](VERIFICATION-RECORD.md)
describes the earlier work and the local package check.

## Reaction table and identifiers

Write concentrations as `(x,y,z) = (x1,x2,x3)` and rates as
`(a,b,c,d,e,f) = (kappa1,...,kappa6)`. The empty complex is `0`.

| Rate | Monomial | Network A | Network B |
|---|---|---|---|
| a | az | X3 → X1 | X3 → X1 + X2 + X3 |
| b | by | X2 → X1 + X3 | X2 → 0 |
| c | cyz | X2 + X3 → 0 | X2 + X3 → X1 |
| d | dx | X1 → 0 | X1 → 0 |
| e | exz | X1 + X3 → X1 + X2 + X3 | X1 + X3 → X1 + X2 |
| f | fxy | X1 + X2 → X1 + X2 + X3 | X1 + X2 → X1 + X2 + X3 |

[reactions.json](reactions.json) is a new transcription of the formulas in
[Example 4, PDF page 22](https://arxiv.org/pdf/2608.01116v1#page=22).
The runner generates `Y` and `N = products - reactants` from these lists; rows
follow species order and columns follow the table. Both matrices `N` have rank 3,
so there are no nonzero linear conservation laws.

| Package label | Example 4 order | Author candidate_index | Author master_index | Position in author list |
|---|---:|---:|---:|---:|
| A | first | 33449 | 106608 | 20 |
| B | second | 239089 | 246746 | 369 |

The identifiers come from the authors' repository at
[commit 323360d](https://github.com/zjdong-sudo/paper-multistability-results/tree/323360df12c0787f9bbf7f223a6229b98233be39).
The list position, called `rank` there, is not matrix rank. No code, JSON dataset
or raw log from that repository is bundled; its redistribution licence was not
established. [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) records the sources.

## Fixed rational samples

The last four rates are `1/2` in every row.

| Sample | a | b | Positive equilibria | Exponentially stable |
|---|---|---|---:|---:|
| A33449-1 | 109/32768 | 235/256 | 3 | 1 |
| A33449-2 | 33/131072 | 31/32 | 3 | 1 |
| B239089-1 | 555/2048 | 103/65536 | 3 | 1 |

[samples.json](samples.json) restates these cited numerical inputs. The two A
vectors were recorded in the authors' logs; B is also printed in Example 4.
The isolation intervals and stability certificates in [certificates/](certificates/)
were computed in this study, not copied from those logs. All nine equilibria are
nondegenerate. Their determinant signs in increasing `z` order are `+, -, +`.
No floating-point approximation decides a root count or stability sign.

The programs are bounded checkers for these fixtures. They are not general solvers
for arbitrary networks, repeated roots or arbitrary parameter vectors. The
all-positive-rate conclusions follow from the two written derivations, rather
than from extrapolating the three examples.

## Verification history and scope

The second derivation and implementation were committed before their author read
the first derivation, implementation or certificate outputs. That check used
`z` as the elimination variable and supplied a direct Routh–Hurwitz argument.
The first used `t=z/y` for A and `y` for B. The comparison matched all nine
equilibria and their stability signs. The task given to the second checker already
stated the expected conclusions and sample counts, so this was not blind testing.
The package preparation used both implementations. No proof-assistant
formalization has been completed.

[provenance.json](provenance.json) gives the research commits and hashes.
The unchanged pre-comparison files are identified separately from the current
runner and its adapted comparison code. Hash checks verify bundled bytes; the
private research commit IDs are provenance references, not public download links
or an independently authenticated timestamp.

This package makes no claim to a new classification, general stability principle
or new elimination method. It specializes Tang–Zhang's Lemma 3 and the general
determinant-and-stability argument of Torres–Feliu, Proposition 2 and Theorem 1.
Novelty or priority of the two-network conclusions has not been established.

Tang–Zhang's Theorem 2 states the upper bound and the 375/373 classification.
The further theoretical explanation of the bound across those networks, and
the common structure of the 373 multistable networks, remain unresolved here.
The package neither proves that full classification nor audits its enumeration.
It does not classify nonlinear asymptotic stability at degenerate equilibria,
boundary attractors, periodic orbits, dissipativity, global dynamics or biological
realizability. The origin is a boundary equilibrium of both networks and is
excluded from every positive-equilibrium count.

## Authorship, sources and acknowledgements

Tatsuhiko Takimoto is responsible for this package. Claude and Codex were used
for derivations, implementation, computation, literature checking and drafting.
Neither the comparison of implementations nor the AI editorial review
constitutes human peer review.

- Xiaoxian Tang and Jiandong Zhang, *The Ubiquity of Three Steady States:
  Minimal Multistable Zero-One Reaction Networks*,
  [arXiv:2608.01116v1](https://arxiv.org/abs/2608.01116v1) (2026), supplied the two
  networks and their exceptional role, the bound of three, the 375/373
  classification and the stability lemma used here. Those are their results.
- Angélica Torres and Elisenda Feliu, *Symbolic proof of bistability in reaction
  networks*, [arXiv:1909.13608v2](https://arxiv.org/abs/1909.13608v2) (2020 version),
  supplied the general reduction-to-determinant and stability-alternation results
  to which these formulas are a specialization.

I am grateful for their work and hope these explicit certificates are useful to
others studying the same networks. No author endorsement is implied.

Code is licensed under Apache-2.0; original text, data and computation records
under CC BY 4.0. See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff).
Cite the original papers for their networks and results.
