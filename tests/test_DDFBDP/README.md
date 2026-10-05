# Specimen DD-FBD tests

Run commands and model scope are documented in `examples/DDFBDP/README.md`.
The integration runner creates an isolated temporary working directory and
preserves console logs alongside the requested JSON report; it does not delete
other test outputs. The morphology example itself writes into output/DDFBDP.

`native_constant.Rev` compares absolute densities (including shape factors and
conditioning) at alpha=0. Fixtures include reversed child order, a terminal
fossil, a fossil ancestor, an artificial ancestor root, an all-fossil tree,
extant-only trees, and decimal fossil ages subject to Newick roundoff. rho is
0.4 or 1 and all applicable origin conditioning options are checked. All-fossil
trees are offset after Newick import to restore their positive youngest age.

`numerical_reference.py` uses Python's standard library and a C++17 compiler.
It compares raw densities before the Rev tree shape factor and conditioning.
Its independent RK4 solver does not call the production propagator. The analytic
reference implements the constant-rate E/q solution. Hidden-count caps 64 and
128 are compared for DD fixtures. It also checks a no-birth ancestor chain and
the impossibility of terminal fossils when mu=0 and rho=1, the complete pure-birth
limit, and finite log densities for extremely rare observations.

Short chains verify integration and proposal/cache behavior; they are not
posterior convergence, parameter recovery, or simulation-based calibration.

Native `dnBDSTP` uses an approximate log-factorial in its shape factor; the new
prior deliberately uses the same label constant. There is one documented
conditioning difference: native `pSurvival(origin,0)` divides the sampled-extant
survival probability by rho. The new prior uses the actual probability. Thus
`lnP_native - lnP_new = ln(rho)` for survival conditioning at alpha=0. The native
regression explicitly corrects this known factor, rather than loosening its
1e-8 tolerance. Time/sampling conditions and rho=1 match directly.

Published console logs replace absolute checkout and temporary-directory paths
with placeholders. Numerical values and diagnostic messages are preserved.
