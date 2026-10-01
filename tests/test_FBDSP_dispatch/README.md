# FBDSP const-virtual dispatch regression

The tree distribution's public `computeLnProbability()` calls a **const** virtual
`computeLnProbabilityDivergenceTimes()`. The former non-const FBDSP method hid that
virtual instead of overriding it, so normal dispatch skipped its range evaluator.

`check_override.cpp` asserts the actual member signature. It fails to compile
against the original header and passes after the fix:

```sh
python3 tests/test_FBDSP_dispatch/check_override.py --boost-include /path/to/boost/include
```

`test_dispatch.cpp` checks the public runtime dispatch through an
`AbstractRootedTreeDistribution&`. It supplies sentinel range/time contributions
in a real FBDSP subclass and checks their sum and the existing tree-shape term.
This isolates dispatch correctness; it does not validate the species-range
mathematics, conditioning, or MCMC behavior.

After building the RevBayes core libraries:

```sh
python3 tests/test_FBDSP_dispatch/check_runtime.py \
  --build-dir /path/to/meson/build \
  --boost-root /path/to/boost/install
```

Both runners honor `CXX` and require C++23. These are explicit compiler/core tests,
not `.Rev` integration scripts discovered by `run_integration_tests.sh`.
