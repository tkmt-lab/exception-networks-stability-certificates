# Verification record

The mathematical inputs are the two networks in Tang–Zhang Example 4 and the
three rational vectors in samples.json. All-parameter proofs and finite-sample
checks are different kinds of evidence.

## Earlier derivations and comparison

The first implementation used sparse integer multivariate polynomials for its
identity checks and `fractions.Fraction` for Sturm sequences, root isolation and
Jacobian interval arithmetic. It used `t=z/y` for A and `y` for B. The files under
`archive/producer/` and its certificates are unchanged copies of that implementation.
Its all-parameter stability argument invokes Tang–Zhang Lemma 3.

The second implementation reconstructed mass-action equations from the input
matrices and eliminated `x,y` with `z` retained. The accompanying derivation
proved positivity of the first two Hurwitz determinants directly from the
equilibrium fluxes. The implementation used SymPy
for identities and a separate Fraction polynomial implementation for Sturm
chains and rational intervals. `check_certificates.py` rechecked the nine roots
using three sign-changing intervals plus cubic degree and 63 zero-free rational
sign tests; it does not import the Fraction polynomial implementation.

The second derivation, code and certificates were fixed at research commit
`9c1dde75aa476ba15bc1bd97c847c2ad2658ad3d` before reading the first derivation or
results. Its task disclosed the proposed conclusions and expected fixture counts.
It shared the network matrices and parameter inputs; it was not independent of
those inputs or blind to the expected result. It checked the reaction diagrams
against the paper. Its author did not read or import the first implementation
during that independent phase or the original result comparison.

The later comparison at `e402dbdf2d2b9484c4aeb2101a15b090b79d8f7f` matched the
symbolic reductions and all nine equilibria. A's root order reverses between
`t` and `z`; B's order agrees between `y` and `z`. The original comparison record
is `certificates/historical-comparison.json`. Its statement that producer code
was unread describes that historical phase, not preparation of this package.

The manifest at `archive/independent/pre-comparison-manifest.json` lists 18 files
from the original research freeze. The package includes only the own-work files
mapped in provenance.json. Source extracts, task instructions and environment
records were omitted. Hashes of omitted files are recorded provenance, not a
claim that their bytes are available or locally checked here. The archived
Japanese derivation is unchanged; docs/z-derivation.md gives its English account.

## What the package command checks

`verify.py` checks the package manifest and the archived-file hashes, reconstructs
the matrices from the new reaction transcription, and executes the unchanged
programs in disposable directories. It compares the recomputed mathematical
fields with the historical certificates. Version, timing and source-file hash
fields that necessarily differ between a new transcription and the historical
input files are not compared as mathematical results.

The final reconciliation adapts the post-freeze comparison code. It verifies
the relation between parametrizations, rational root counts, the mapping of
root intervals to coordinates, and matching Hurwitz signs. This is a rerun
using both implementations, not a third independent derivation.

Root-count and sign decisions use exact rational arithmetic. Floating point is
used only for elapsed-time measurement. No assertion is made that the Python
interpreter, SymPy or either checker has itself been formally verified. Hashes
detect changed bytes relative to the bundled manifest; they do not authenticate
the historical date or independently prove the manifest's trustworthiness.

No human expert review or proof-assistant formalization has been completed.
The prose arguments for complete positive branches, inequalities and the
all-parameter root bound remain part of the mathematical evidence to inspect.

## Local package check

On 2026-10-24 JST, a local clone of package commit
`81e9df591df0fbfea611580893f748dcacb885fa` passed the command below in a newly
created virtual environment. This is the initial package commit, before the
editorial clarifications to README.md and this record.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install --disable-pip-version-check -r requirements.txt
.venv/bin/python -B verify.py
```

The environment was macOS, Python 3.14.7, SymPy 1.14.0 and mpmath 1.3.0.
Verification took 20.83 seconds and returned exit status 0. The complete clone,
installation and verification took 26.49 seconds, with maximum resident set size
173,998,080 bytes. These are maintainer-reported measurements, not a separate
review. The original research workspace supplied no runtime input files.

The final preparation report, including the subsequent check of the final
package commit, is kept in the research workspace and is not bundled here.
GitHub CI is configured to run the same verification command on Ubuntu; it has
not been run remotely at this preparation stage.
