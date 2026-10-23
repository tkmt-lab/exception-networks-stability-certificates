"""Exact integer polynomial checks for the two R023 reductions; no CAS dependency.

This is the producer's algebra check, NOT an independent verifier or root certificate.
Run: python3 -B studies/R023/check_identities.py
"""
from collections import defaultdict
from pathlib import Path
import hashlib
import itertools
import json
import resource
import time

NAMES = "a b c d e f t x y z".split()
NV = len(NAMES)


class P:
    def __init__(self, value=0):
        if isinstance(value, P):
            self.terms = value.terms
        elif isinstance(value, dict):
            self.terms = {k: v for k, v in value.items() if v}
        else:
            self.terms = {(0,) * NV: value} if value else {}

    def __add__(self, other):
        out = defaultdict(int, self.terms)
        for k, v in P(other).terms.items():
            out[k] += v
        return P(out)

    __radd__ = __add__

    def __neg__(self):
        return P({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + -P(other)

    def __rsub__(self, other):
        return P(other) + -self

    def __mul__(self, other):
        out = defaultdict(int)
        for ka, va in self.terms.items():
            for kb, vb in P(other).terms.items():
                out[tuple(a + b for a, b in zip(ka, kb))] += va * vb
        return P(out)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        assert isinstance(exponent, int) and exponent >= 0
        out = P(1)
        for _ in range(exponent):
            out = out * self
        return out

    def derivative(self, name):
        i = NAMES.index(name)
        out = {}
        for k, v in self.terms.items():
            if k[i]:
                key = list(k)
                key[i] -= 1
                out[tuple(key)] = v * k[i]
        return P(out)

    def clear_substitution(self, numerators, denominator, degree):
        """Return D**degree * P(a,..,f,t,X/D,Y/D,Z/D), exactly."""
        out = P()
        for k, v in self.terms.items():
            spatial_degree = sum(k[7:])
            assert spatial_degree <= degree
            term = P({k[:7] + (0, 0, 0): v})
            for numerator, power in zip(numerators, k[7:]):
                term *= numerator ** power
            out += term * denominator ** (degree - spatial_degree)
        return out

    def digest(self):
        serialized = json.dumps(sorted(self.terms.items()), separators=(",", ":"))
        return hashlib.sha256(serialized.encode()).hexdigest()


def variable(name):
    k = [0] * NV
    k[NAMES.index(name)] = 1
    return P({tuple(k): 1})


def determinant(m):
    return (m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1])
            - m[0][1] * (m[1][0]*m[2][2] - m[1][2]*m[2][0])
            + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0]))


def run():
    start = time.monotonic()
    here = Path(__file__).resolve().parent
    spec = json.loads((here / "networks.json").read_text())
    a, b, c, d, e, f, t, x, y, z = [variable(n) for n in NAMES]
    variables = [x, y, z]
    rates = [a, b, c, d, e, f]
    results = []
    for row in spec["networks"]:
        ident, Y, N = row["candidate_index"], row["Y"], row["N"]
        ode = [P() for _ in range(3)]
        for j in range(6):
            rate = rates[j]
            for i in range(3):
                rate *= variables[i] ** Y[i][j]
            for i in range(3):
                ode[i] += N[i][j] * rate
        if ident == 33449:
            expected = [-d*x+b*y+a*z, -b*y-c*y*z+e*x*z, -a*z+b*y-c*y*z+f*x*y]
            L = e*a*t + e*b - c*d
            D = t*L
            nums = [b*(b+a*t), b*d, b*d*t]
            q = -e*a*a*t**3 + a*c*d*t**2 + b*(e*b+f*a-2*c*d)*t + b*b*f
            residuals = [P(), P(), b*d*q]
            det_expected = -b*d*D*(q.derivative("t")*D - 2*q*D.derivative("t"))
            branch = "t>max(0,(c*d-e*b)/(e*a)); D=t*(e*a*t+e*b-c*d)>0"
        elif ident == 239089:
            expected = [-d*x+a*z+c*y*z, -b*y+a*z+e*x*z-c*y*z, f*x*y-e*x*z-c*y*z]
            A = a+c*y
            U = f*A-c*d
            D = d*e*A
            nums = [y*U*A, y*D, d*y*U]
            q = -b*d*e*A + d*(a-c*y)*U + y*U**2
            expanded = c*c*f*f*y**3 + c*f*(2*a*f-3*c*d)*y**2 + ((a*f-c*d)**2+c*c*d*d-b*c*d*e)*y + a*d*(a*f-c*d-b*e)
            assert not (q-expanded).terms
            residuals = [P(), y*q*D, P()]
            det_expected = y*U*D*((q+y*q.derivative("y"))*D-y*q*D.derivative("y"))
            branch = "y>max(0,(c*d-a*f)/(c*f)); U=f*(a+c*y)-c*d>0"
        else:
            raise ValueError(ident)
        checks = {}
        for i in range(3):
            checks[f"ode_{i+1}"] = not (ode[i]-expected[i]).terms
            checks[f"substitution_{i+1}"] = not (ode[i].clear_substitution(nums, D, 2)-residuals[i]).terms
        jacobian = [[r.derivative(n) for n in ("x", "y", "z")] for r in ode]
        scaled_j = [[p.clear_substitution(nums, D, 1) for p in r] for r in jacobian]
        checks["determinant_identity"] = not (determinant(scaled_j)-det_expected).terms
        minors = []
        for cols in itertools.combinations(range(6), 3):
            det = determinant([[N[i][j] for j in cols] for i in range(3)])
            if det:
                minors.append({"columns_one_based": [j+1 for j in cols], "determinant": det})
        checks["rank_three"] = bool(minors)
        assert all(checks.values()), (ident, checks)
        results.append({"candidate_index": ident, "checks": checks, "rank_witness": minors[0],
                        "q_polynomial_sha256": q.digest(), "positive_branch": branch})
    return {"status": "PASS", "arithmetic": "exact integer coefficients, sparse multivariate polynomials",
            "limitations": "Same-agent producer check; no independent review, no root isolation or global stability proof checked here.",
            "results": results, "elapsed_seconds": round(time.monotonic()-start, 6),
            "max_rss_platform_units": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "rss_unit_macos": "bytes", "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "networks_sha256": hashlib.sha256((here / "networks.json").read_bytes()).hexdigest()}


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, indent=2))
