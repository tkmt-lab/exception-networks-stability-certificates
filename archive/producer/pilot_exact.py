"""Bounded exact pilot: three fixed published samples, rational Sturm/RH only.

This is a same-agent producer prototype, not the independent checker in plan.md.
No external dependencies. Stops on rational split roots or unresolved signs.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import resource
import time


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def derivative(p):
    return trim([i*p[i] for i in range(1, len(p))] or [F(0)])


def remainder(p, q):
    p = p[:]
    while len(p) >= len(q) and p != [0]:
        shift, factor = len(p)-len(q), p[-1]/q[-1]
        for i, value in enumerate(q):
            p[i+shift] -= factor*value
        trim(p)
    return p


def evaluate(p, x):
    out = F(0)
    for coefficient in reversed(p):
        out = out*x + coefficient
    return out


def sturm(p):
    out = [p, derivative(p)]
    while len(out[-1]) > 1:
        r = [-v for v in remainder(out[-2], out[-1])]
        assert r != [0], "Repeated roots: handle separately"
        out.append(r)
    return out


def variations(seq, x):
    signs = [1 if v > 0 else -1 for p in seq if (v := evaluate(p, x)) != 0]
    return sum(a != b for a, b in zip(signs, signs[1:]))


class I:
    def __init__(self, lo, hi=None):
        if isinstance(lo, I):
            self.lo, self.hi = lo.lo, lo.hi
        else:
            self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, b):
        b = I(b)
        return I(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, b):
        return self + -I(b)

    def __rsub__(self, b):
        return I(b) + -self

    def __mul__(self, b):
        b = I(b)
        vals = [a*c for a in (self.lo, self.hi) for c in (b.lo, b.hi)]
        return I(min(vals), max(vals))

    __rmul__ = __mul__

    def __truediv__(self, b):
        b = I(b)
        assert b.lo > 0 or b.hi < 0, "Denominator interval contains zero"
        return self * I(1/b.hi, 1/b.lo)

    def json(self):
        return [str(self.lo), str(self.hi)]


def det(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
            - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))


def jacobian(row, rates, coordinates):
    Y, N = row["Y"], row["N"]
    J = [[I(0) for _ in range(3)] for _ in range(3)]
    for j in range(6):
        for variable in range(3):
            if Y[variable][j] == 0:
                continue
            drate = I(rates[j])
            for species in range(3):
                power = Y[species][j] - (species == variable)
                for _ in range(power):
                    drate *= coordinates[species]
            for equation in range(3):
                J[equation][variable] += N[equation][j] * drate
    return J


def run():
    start = time.monotonic()
    deadline = start + 60
    here = Path(__file__).resolve().parent
    spec = json.loads((here / "networks.json").read_text())
    samples = json.loads((here / "samples.json").read_text())
    rows = {r["candidate_index"]: r for r in spec["networks"]}
    results = []
    for sample in samples["samples"]:
        ident = sample["candidate_index"]
        rates = [F(s) for s in sample["kappa"]]
        a, b, c, d, e, f = rates
        if ident == 33449:
            q = [b*b*f, b*(e*b+f*a-2*c*d), a*c*d, -e*a*a]
            lower = max(F(0), (c*d-e*b)/(e*a))

            def lift(t):
                L = e*a*t+e*b-c*d
                D = t*L
                return [I(b)*(b+a*t)/D, I(b*d)/D, I(b*d)/L]
        else:
            q = [a*d*(a*f-c*d-b*e), (a*f-c*d)**2+c*c*d*d-b*c*d*e,
                 c*f*(2*a*f-3*c*d), c*c*f*f]
            lower = max(F(0), (c*d-a*f)/(c*f))

            def lift(y):
                A = a+c*y
                U = f*A-c*d
                return [y*U/(d*e), y, y*U/(e*A)]

        seq = sturm(q)
        # Strict Cauchy bound: every root has |root| < upper.
        upper = 2 + max(abs(coef/q[-1]) for coef in q[:-1])
        assert lower < upper and evaluate(q, lower) != 0 and evaluate(q, upper) != 0
        count = lambda lo, hi: variations(seq, lo)-variations(seq, hi)
        total = count(lower, upper)
        pending = [(lower, upper)]
        isolated = []
        while pending:
            assert time.monotonic() < deadline, "60-second pilot cap"
            lo, hi = pending.pop()
            n = count(lo, hi)
            if n == 0:
                continue
            if n == 1 and hi-lo < F(1, 2**48):
                isolated.append((lo, hi))
                continue
            mid = (lo+hi)/2
            assert evaluate(q, mid) != 0, "Rational split root: separate branch required"
            pending.extend([(lo, mid), (mid, hi)])
        roots = []
        for lo, hi in sorted(isolated):
            coords = lift(I(lo, hi))
            assert all(c.lo > 0 for c in coords)
            J = jacobian(rows[ident], rates, coords)
            b1 = -sum(J[i][i] for i in range(3))
            b2 = sum(J[i][i]*J[j][j]-J[i][j]*J[j][i] for i, j in [(0, 1), (0, 2), (1, 2)])
            b3 = -det(J)
            H2 = b1*b2-b3
            if all(v.lo > 0 for v in [b1, b2, b3, H2]):
                verdict = "stable_Routh_Hurwitz"
            elif b3.hi < 0:
                verdict = "unstable_positive_real_eigenvalue"
            else:
                raise AssertionError("Unresolved rational interval signs")
            roots.append({"interval": [str(lo), str(hi)], "variations": [variations(seq,lo),variations(seq,hi)],
                          "coordinates": [c.json() for c in coords],
                          "characteristic_b1_b2_b3_H2": [v.json() for v in [b1,b2,b3,H2]], "verdict": verdict})
        assert len(roots) == total
        results.append({"sample": sample["id"], "candidate_index": ident, "q_ascending": list(map(str,q)),
                        "admissible_lower": str(lower), "strict_cauchy_upper": str(upper),
                        "total_variations": [variations(seq,lower),variations(seq,upper)],
                        "distinct_positive_equilibria": total, "squarefree": len(seq[-1]) == 1,
                        "stable": sum(r["verdict"] == "stable_Routh_Hurwitz" for r in roots), "roots": roots})
    return {"status": "PASS", "arithmetic": "fractions.Fraction only for every mathematical decision",
            "limitations": "Producer pilot, depends on the separately checked reduction. No independent checker or external review yet.",
            "elapsed_seconds": round(time.monotonic()-start,6),
            "max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "inputs_sha256": {p: hashlib.sha256((here/p).read_bytes()).hexdigest() for p in ["networks.json","samples.json"]},
            "results": results}


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, indent=2))
