# Kernel regression tests

`pytest tests/regression/` runs the current GraKeL kernels against fixed output
matrices. It detects unintended behaviour changes; it does not install or run
historical GraKeL releases.

    data/graphs.json              fixed corpus
    data/specifications.json      graph split and seven kernel configurations
    expected/                     frozen train/test matrices and tolerances
    test_regression.py            the only test

The references were captured from GraKeL 0.1.10. Most comparisons are exact to
`1e-9`; NSPD uses `rtol=0.05`, `atol=0.0001` because its implementation changed
again after 0.1.10. Update a frozen matrix only after intentionally reviewing a
kernel behaviour change.
