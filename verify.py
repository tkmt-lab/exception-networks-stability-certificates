#!/usr/bin/env python3
"""Check file integrity, rerun both derivations, and compare fixed certificates.

SPDX-License-Identifier: Apache-2.0
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def same(actual, expected, label):
    # JSON serialization distinguishes true from 1 as well as checking values.
    require(json.dumps(actual, sort_keys=True) == json.dumps(expected, sort_keys=True),
            f"fixed mathematical certificate differs: {label}")


def integrity():
    manifest = read(ROOT / "manifest.json")
    for name, digest in manifest["files"].items():
        path = ROOT / name
        require(path.is_file() and not path.is_symlink(), f"missing file: {name}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == digest,
                f"file hash mismatch: {name}")
    for item in read(ROOT / "provenance.json")["historical_files"]:
        require(manifest["files"][item["package_file"]] == item["sha256"],
                f"historical file changed: {item['package_file']}")
    claims = read(ROOT / "claims.json")
    require(claims["prepared_not_released"] is False, "release-state metadata differs")
    require(claims["network_ids"] == [33449, 239089], "claim network set differs")
    require(claims["fixed_samples"] == [
        {"id": name, "positive_equilibria": 3, "exponentially_stable": 1}
        for name in ["A33449-1", "A33449-2", "B239089-1"]], "fixture claims differ")
    zenodo = read(ROOT / ".zenodo.json")
    citation = read(ROOT / "CITATION.cff")  # JSON is a YAML 1.2 subset.
    require(zenodo["upload_type"] == citation["type"] == "dataset", "dataset metadata differs")
    require(zenodo["title"] == citation["title"] == claims["title"], "titles differ")
    require(zenodo["creators"][0]["name"] == "Tatsuhiko Takimoto", "creator differs")


def inputs():
    data = read(ROOT / "reactions.json")
    species = data["species"]
    require(species == ["X1", "X2", "X3"], "species order differs")
    require([n["candidate_index"] for n in data["networks"]] == [33449, 239089],
            "network case set differs")
    networks = []
    for network in data["networks"]:
        reactions = network["reactions"]
        require(len(reactions) == 6, "six ordered reactions are required")
        for reaction in reactions:
            require(set(reaction) == {"reactants", "products"}, "invalid reaction record")
            for side in reaction.values():
                require(len(set(side)) == len(side) and set(side) <= set(species),
                        "a zero-one species list is required")
        Y = [[r["reactants"].count(s) for r in reactions] for s in species]
        N = [[r["products"].count(s) - r["reactants"].count(s) for r in reactions]
             for s in species]
        networks.append({"candidate_index": network["candidate_index"], "Y": Y, "N": N})
    samples = read(ROOT / "samples.json")["samples"]
    expected = read(ROOT / "certificates/independent-samples.json")["samples"]
    # Source descriptions may change; the three parameter vectors may not.
    fields = ("id", "candidate_index", "kappa")
    sample_input = [{key: sample[key] for key in fields} for sample in samples]
    same(sample_input, [{key: sample[key] for key in fields} for sample in expected], "sample inputs")
    require(all(len(s["kappa"]) == 6 and all(Fraction(k) > 0 for k in s["kappa"])
                for s in samples), "strictly positive rational rates required")
    return {"networks": networks}, {"samples": sample_input}


def execute(stage, script):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    # Assertions in the unchanged historical scripts must remain enabled.
    env.pop("PYTHONOPTIMIZE", None)
    env.pop("PYTHONPATH", None)
    result = subprocess.run([sys.executable, "-E", "-B", str(stage / script)],
                            cwd=stage, env=env, capture_output=True, text=True, timeout=120)
    require(result.returncode == 0, f"{script} failed:\n{result.stdout}\n{result.stderr}")
    return result.stdout


def main():
    start = time.monotonic()
    require(not sys.flags.optimize, "run without Python -O; historical checkers use assertions")
    import sympy
    import mpmath
    require(sympy.__version__ == "1.14.0" and mpmath.__version__ == "1.3.0",
            "install the versions in requirements.txt")
    integrity()
    networks, samples = inputs()
    # All generated files live in a disposable copy; saved certificates are read-only inputs.
    with tempfile.TemporaryDirectory(prefix=".verify-", dir=ROOT) as temporary:
        work = Path(temporary)
        producer, independent = work / "producer", work / "independent"
        shutil.copytree(ROOT / "archive/producer", producer)
        shutil.copytree(ROOT / "archive/independent", independent)
        (independent / "input").mkdir()
        for directory in [producer, independent / "input"]:
            (directory / "networks.json").write_text(json.dumps(networks, indent=2) + "\n")
            (directory / "samples.json").write_text(json.dumps(samples, indent=2) + "\n")
        for script, certificate in [("check_identities.py", "identity-check.json"),
                                    ("pilot_exact.py", "pilot-exact.json")]:
            result = json.loads(execute(producer, script))
            same(result["results"], read(ROOT / "certificates" / certificate)["results"], certificate)
            (independent / certificate).write_text(json.dumps(result) + "\n")
        print("PASS: first derivation and three exact rational fixtures", flush=True)
        for script, certificate, key in [
            ("derive_checker.py", "independent-identities.json", "networks"),
            ("certify_samples.py", "independent-samples.json", "samples"),
            ("check_certificates.py", "certificate-recheck.json", "samples")]:
            execute(independent, script)
            same(read(independent / certificate)[key],
                 read(ROOT / "certificates" / certificate)[key], certificate)
        print("PASS: second derivation, Fraction Sturm certificates and separate sign checks", flush=True)
        shutil.copyfile(ROOT / "tools/reconcile.py", independent / "reconcile.py")
        execute(independent, "reconcile.py")
        comparison = read(independent / "comparison-check.json")
        historical = read(ROOT / "certificates/historical-comparison.json")
        for key in ["status", "symbolic_comparison", "samples", "mathematical_discrepancies"]:
            same(comparison[key], historical[key], "reconciliation: " + key)
        # Tie the advertised summary to recomputed values, not only to source hashes.
        for claim, actual in zip(read(ROOT / "claims.json")["fixed_samples"],
                                 read(independent / "independent-samples.json")["samples"]):
            same(claim, {"id": actual["id"], "positive_equilibria": actual["positive_equilibria"],
                         "exponentially_stable": actual["exponentially_stable_equilibria"]}, "claim counts")
    print("PASS: 2 networks, 3 samples, 9 equilibria; both implementations agree")
    print(f"Elapsed: {time.monotonic() - start:.2f} seconds")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, ValueError, KeyError, OSError, ImportError,
            subprocess.TimeoutExpired) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
