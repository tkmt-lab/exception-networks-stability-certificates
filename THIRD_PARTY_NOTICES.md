# Sources and included material

The reactions in `reactions.json` and README.md were newly transcribed from
Tang–Zhang, *The Ubiquity of Three Steady States: Minimal Multistable Zero-One
Reaction Networks*, [arXiv:2608.01116v1, Example 4, page 22](https://arxiv.org/pdf/2608.01116v1#page=22).
The mathematical reaction formulas are restated in a species-list representation;
the original figure, paper text and PDF are not reproduced.

The candidate/master identifiers and list positions were matched against the
authors' [repository at commit 323360df12c0787f9bbf7f223a6229b98233be39](https://github.com/zjdong-sudo/paper-multistability-results/tree/323360df12c0787f9bbf7f223a6229b98233be39),
in `results/375_multistable_networks.json`. Its `rank` field is list position.
No file from that repository is redistributed.

The first two parameter vectors in `samples.json` restate the rational numbers
reported in `data/check_multistability_2.txt`, lines 399–402 and 414–415, at that
commit. The third restates the vector printed in Example 4 and also found in
`data/check_multistability_4.txt`, lines 1003–1006. Direct links appear beside
each vector. These numerical facts are inputs, not author-repository log extracts.
The original logs, computed author outputs, data tables and source code are absent.
No redistribution licence for that repository was established in the source check.

Every bundled Python checker and computed certificate was written or produced in
the R023 study. The Apache-2.0 and CC BY 4.0 notices cover that original work;
they do not grant rights in linked third-party material or claim ownership of
the mathematical facts being restated.

SymPy 1.14.0 and mpmath 1.3.0 are separately installed dependencies, not vendored
code. Their own licences apply to their distributions. The full Apache-2.0 and
CC BY 4.0 texts are included for the package's licensing split.
