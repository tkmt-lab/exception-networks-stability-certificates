# Derivation with z retained

This is an English account of the independently written derivation fixed at
research commit `9c1dde75aa476ba15bc1bd97c847c2ad2658ad3d`. Its unchanged Japanese
text is [archived here](../archive/independent/derivation.md). The expected
conclusions and sample counts were known to the checker; the first derivation,
code and results had not been read. All six rates are strictly positive.

Use the mass-action equations and rate order in [the first derivation](ratio-y-derivation.md).
They can also be reconstructed independently from `N diag(a,...,f)` times
the monomials `(z,y,yz,x,xz,xy)`. Both stoichiometric ranks are 3.

## Complete parametrizations of the positive equilibria

For A, solving `F1=F2=0` as linear equations for `x,y` gives

```text
D_A(z) = bd+(cd-be)z,
x(z) = az(b+cz)/D_A(z),
y(z) = ae z^2/D_A(z),
I_A = {z>0 : D_A(z)>0}.
```

The interval is `(0,infinity)` if `cd>=be`, and `(0,bd/(be-cd))` otherwise.
At a positive equilibrium, the same linear equations imply
`D_A x=az(b+cz)>0`, so neither `D_A=0` nor `D_A<0` can hide a positive solution.

Put

```text
P_A(z) = D_A(z)^2 - ez(b-cz)D_A(z) - aef z^2(b+cz).
```

The remaining equation is `F3(x(z),y(z),z)=-az P_A/D_A^2`. Thus roots of
`P_A` in `I_A` are in bijection with positive equilibria. The degree is at
most three, the constant term is `b^2 d^2>0`, and the cubic coefficient is
`ce(cd-be-af)`. That coefficient may vanish at positive rates. For example
`(a,b,c,d,e,f)=(1,1,1,2,1,1)` gives `P_A=z^2+2z+4`.

For B the corresponding formulas are

```text
D_B(z) = d(b+cz)-ce z^2,
x(z) = az(b+2cz)/D_B(z),
y(z) = az(d+ez)/D_B(z),
I_B = {z>0 : D_B(z)>0},
P_B(z) = D_B(z)(cd+eb+3cez)-af(b+2cz)(d+ez).
```

Here `I_B=(0,rho)`, where
`rho=(cd+sqrt(c^2 d^2+4bcde))/(2ce)`. At a positive equilibrium,
`D_B x=az(b+2cz)>0`, so the omitted denominator branches cannot be positive.
On this curve `F3=-az^2 P_B/D_B^2`. Again the roots in the indicated interval
are exactly the positive equilibria. The leading coefficient is
`-3c^2 e^2<0`, so B's polynomial is always cubic. Its constant term is
`bd(cd+eb-af)`; a possible root at `z=0` is not a positive equilibrium.

## Determinant identities

For both models, the minor of `(F1,F2)` with respect to `(x,y)` equals `D`.
On the parametrized curve the chain rule and a Schur complement give

```text
det J = D d/dz [F3(x(z),y(z),z)].
```

At admissible roots this becomes

```text
det J_A = -az P_A'(z)/D_A(z),
det J_B = -az^2 P_B'(z)/D_B(z).
```

The prefactors are strictly negative. In particular degeneracy is equivalent
to a repeated root. The symbolic checker verifies the full off-root chain-rule
identity, as well as its root-restricted consequence.

This determinant reduction is a direct instance of Torres–Feliu,
[Proposition 2](https://arxiv.org/abs/1909.13608v2): here the stoichiometric rank
and number of species are both 3, there are no conservation rows, the retained
equation and coordinate are both the third, and the minor is `D>0`.
Their Theorem 1 uses the Hurwitz positivity established below. Completeness of
the positive parametrizations is supplied by the preceding argument, not by
the determinant identity alone.

## Direct Routh–Hurwitz argument

Let `v` be the six positive equilibrium reaction fluxes. Solving `Nv=0` gives

```text
A: v = (p,q,r,p+q,q+r,p-q+r),
B: v = (p,p-q+r,q,p+q,r,q+r),
```

where `p,q,r>0` and `p-q+r>0`. The flux Jacobians are

```text
J_A = [ -(p+q)/x      q/y          p/z      ]
      [  (q+r)/x    -(q+r)/y       q/z      ]
      [ (p-q+r)/x     p/y        -(p+r)/z   ]

J_B = [ -(p+q)/x      q/y         (p+q)/z   ]
      [     r/x     -(p+r)/y   (p-q+r)/z   ]
      [     q/x        r/y       -(q+r)/z  ].
```

All diagonal entries are negative and all off-diagonal entries are positive.
The principal minors of order two are strictly positive:

| Network | M12 | M13 | M23 |
|---|---|---|---|
| A | p(q+r)/(xy) | q(2p+r)/(xz) | r(p+q+r)/(yz) |
| B | p(p+q+r)/(xy) | r(p+q)/(xz) | q(p+2r)/(yz) |

Write the diagonal as `-A0,-B0,-C0`, and the characteristic polynomial as
`lambda^3+b1 lambda^2+b2 lambda+b3`. Then `b1=A0+B0+C0>0`,
`b2=M12+M13+M23>0`, and direct expansion gives

```text
H2 = b1 b2 - b3
   = (A0+B0)M12 + (A0+C0)M13 + (B0+C0)M23
     + 2A0 B0 C0 + J12 J23 J31 + J13 J21 J32 > 0.
```

The cubic Routh–Hurwitz conditions therefore reduce to `b3=-det J>0`.
A negative determinant gives exponential stability. A positive determinant
makes `b3<0`, so the characteristic polynomial has a positive real root and
the equilibrium is unstable. A zero determinant gives a zero eigenvalue,
which excludes exponential stability; nonlinear asymptotic stability is not
classified. This proves directly the specialization of Tang–Zhang Lemma 3
needed for these two networks.

## Consequences for every positive rate choice

The nonzero polynomials have degree at most three and the parametrizations
are bijections, giving at most three distinct positive equilibria.

If three exist, the polynomials have three distinct roots `z1<z2<z3` in the
admissible interval. They are then cubic and all roots are simple. B has negative
leading coefficient. For A, writing `P_A=L(z-z1)(z-z2)(z-z3)` and using
`P_A(0)>0` also forces `L<0`. The derivative signs are therefore `-,+,-`,
and the determinant signs are `+,-,+`. Exactly the middle equilibrium is
exponentially stable.

More generally, exponential stability requires `P'>0` at an admissible root.
For A, suppose two positive roots `u<v` had positive derivative. The polynomial
would be negative just before each and positive just after each. Since
`P_A(0)>0`, the intermediate value theorem forces one further root in `(0,u)`
and another in `(u,v)`. Four distinct roots contradict degree at most three.
For B, `P_B'` is a quadratic with negative leading coefficient, so its positive
set is one interval or is empty. The polynomial is strictly increasing there
and has at most one zero. Thus each network has at most one exponentially
stable positive equilibrium even when other positive equilibria are degenerate.

## Boundary and verification limits

The only equilibrium on the boundary of the nonnegative orthant is the origin.
For A, `x=0` in `F1` forces `y=z=0`; `y=0` in `F2` forces `xz=0`; and
`z=0` in `F2` forces `y=0`, with the remaining equations giving the origin.
For B, `x=0` in `F1` forces `z=0`; `y=0` in `F2` forces `z=0`; and
`z=0` in `F1,F2` forces `x=y=0`.

At the origin A's characteristic polynomial is `(lambda+a)(lambda+b)(lambda+d)`;
B's is `lambda(lambda+b)(lambda+d)`. Boundary dynamics are not included in the
positive-equilibrium count. In particular the absence-of-boundary-equilibria
assumption in Conradi–Feliu–Mincheva–Wiuf,
[Theorem 1](https://arxiv.org/abs/1608.03993v2), fails, so that theorem is not
applied here. This does not prevent use of the Torres–Feliu reduction above.

The symbolic program checks the ODEs, elimination, residuals, determinant chain
rule, flux Jacobians, principal minors and Hurwitz identity. The branch and
root-sign arguments are prose proofs. Neither this derivation nor the package
has been formally verified or reviewed by a human specialist.
