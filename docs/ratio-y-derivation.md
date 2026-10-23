# Derivation using t=z/y for A and y for B

This is an English presentation of the first R023 derivation (`setup.md` at the
input research commit recorded in provenance.json). It concerns only the two
networks in [Tang–Zhang, Example 4](https://arxiv.org/pdf/2608.01116v1#page=22).
All rates `a,b,c,d,e,f` are strictly positive, in the order of README.md.

The mass-action equations obtained by summing each reaction's stoichiometric
change times its rate are

```text
A: F1 = az + by - dx
   F2 = -by - cyz + exz
   F3 = -az + by - cyz + fxy

B: F1 = az + cyz - dx
   F2 = az - by - cyz + exz
   F3 = -cyz - exz + fxy.
```

The first three columns of the stoichiometric matrices have determinants 3
and 1, respectively. Thus both ranks are 3 and no conservation-law reduction
of the Jacobian is needed. The polynomial checker reconstructs these equations
from the matrices generated from the reaction lists.

## A: the ratio t=z/y

At a positive equilibrium put `t=z/y>0`. Define

```text
L(t) = ea t + eb - cd,
D(t) = t L(t),
T_A = max(0, (cd-eb)/(ea)).
```

From `F1=0`, `x=y(b+at)/d`. Substituting into `F2=0` and dividing by `y`
gives `y t L=bd`. Therefore `L>0`, and every positive equilibrium has

```text
t > T_A,
x = b(b+at)/D,   y = bd/D,   z = bd/L.
```

Conversely these coordinates are positive and satisfy `F1=F2=0` throughout
that interval. There is no missing positive branch at `L=0`, since `y t L=bd>0`.
The excluded cases `y=0` or `t=0` are outside the positive orthant.

With

```text
Q_A(t) = -ea^2 t^3 + acd t^2 + b(eb+fa-2cd)t + b^2 f,
```

substitution gives `F3=bd Q_A/D^2`. Positive equilibria and roots of `Q_A`
in `(T_A,infinity)` are consequently in bijection. Its leading coefficient
is strictly negative for all the allowed rates.

For the determinant, change input coordinates from `(x,y,z)` to `(x,y,t)`
with `z=ty`; the coordinate-change determinant is `y`. At the parametrized
points, the minor of `(F1,F2)` with respect to `(x,y)` is `-yD`.
The chain rule along `F1=F2=0` then gives the off-root identity

```text
D^3 det J_A = -bd D (Q_A' D - 2 Q_A D').
```

At a root this reduces to `det J_A = -y Q_A'(t)`.

## B: the coordinate y

Define

```text
A(y) = a+cy,
U(y) = f(a+cy)-cd,
T_B = max(0, (cd-af)/(cf)),
D(y) = de(a+cy).
```

From `F1=0`, `x=A z/d`. Substituting into `F3=0` and dividing by `z>0`
gives `z=yU/(eA)` and `x=yU/(de)`. Positivity is equivalent to `y>T_B`.
All denominators are positive there. The branch `U=0` has `x=z=0` and cannot
contain a positive equilibrium.

```text
Q_B(y) = -bde(a+cy) + d(a-cy)U + yU^2
       = c^2 f^2 y^3 + cf(2af-3cd)y^2
         + [(af-cd)^2+c^2 d^2-bcde]y + ad(af-cd-be).
```

The parametrization satisfies `F1=F3=0` and `F2=y Q_B/D`. It gives a
bijection between positive equilibria and roots of `Q_B` in `(T_B,infinity)`.
This cubic has strictly positive leading coefficient.

The minor of `(F1,F3)` with respect to `(x,z)` is `yU`. Reordering both inputs
and outputs introduces two sign reversals, which cancel. Differentiating the
remaining equation gives

```text
D^3 det J_B = yU D [(Q_B+y Q_B')D - y Q_B D'].
```

At a root, `det J_B = y^2 U Q_B'(y)/D`.

## Counts and stability

Both nonzero polynomials have degree exactly three for all strictly positive
rates. Their admissible roots therefore give at most three distinct positive
equilibria. Nondegeneracy is equivalent to simplicity of the corresponding
root, because the determinant prefactors are nonzero on the admissible domains.

If three distinct positive equilibria exist, they exhaust the cubic's roots,
so all three are simple. In increasing reduced-coordinate order the derivatives
of `Q_A` have signs `-,+,-`, and those of `Q_B` have signs `+,-,+`. The
Jacobian determinants have signs `+,-,+` for both networks.

[Tang–Zhang Lemma 3](https://arxiv.org/pdf/2608.01116v1#page=8) gives exponential
stability exactly when `det J<0` at a nondegenerate positive equilibrium of
these zero-one networks. Hence exactly the middle point is exponentially stable.
The separate [z derivation](z-derivation.md) proves the required Hurwitz
positivity directly for these two networks as well.

For A, exponential stability requires `Q_A'>0`. The derivative is a quadratic
with negative leading coefficient; its positive set is an interval or is empty.
`Q_A` is strictly increasing there and has at most one zero. For B, stability
requires `Q_B'<0`; the derivative has positive leading coefficient and its
negative set is likewise an interval or is empty. Again there is at most one
zero there. A repeated root gives a zero eigenvalue and cannot be exponentially
stable. This proves the bound of one even when degenerate equilibria coexist.

The argument does not decide nonlinear asymptotic stability at repeated roots.
It does not assert that every rate choice has three positive equilibria.

## Relation to earlier theory and checks

The determinant identity after solving all but one equilibrium equation is
the mechanism of Torres–Feliu, [Proposition 2](https://arxiv.org/abs/1909.13608v2).
Their Theorem 1 derives alternating stability from the determinant identity
and Hurwitz positivity. Here the additional network-specific work is the
explicit parametrization, its positive domain and the cubic coefficient signs.
This is a specialization of existing results, not a new elimination method.

`archive/producer/check_identities.py` checks the ODEs, substitutions, ranks and
the two off-root determinant identities using sparse integer polynomials.
`archive/producer/pilot_exact.py` applies exact rational Sturm and interval
arithmetic to the three fixed vectors. Branch completeness and the root-sign
arguments above are written mathematics, not machine-checked formal proofs.
