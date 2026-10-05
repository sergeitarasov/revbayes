# Specimen diversity-dependent fossilized birth-death prior

This branch implements **`dnDDFBDP`**, also available as
`dnDiversityDependentSpecimenFBD`. It is a prior on an ordinary sampled-ancestor
`TimeTree`, which is supplied directly to `dnPhyloCTMC` for joint inference of the
tree and model parameters. The runnable example uses morphology only.

## Model

The per-lineage speciation rate is

```text
lambda(N) = lambda0 * exp(-alpha * (N - 1)), alpha >= 0
```

`N` counts all living lineages, including unobserved lineages. `lambda0` is the
rate at diversity one; larger `alpha` means stronger decline. Extinction `mu`
and fossil preservation/sampling `psi` are constant per-lineage rates. The total
birth hazard is `N * lambda(N)`. There is no hard carrying capacity or skyline.

Each taxon label represents one fossil occurrence or one living observation.
A fossil can be a terminal tip or a sampled ancestor. Sampling does not kill the
lineage. After a terminal fossil, its continuation becomes hidden rather than
being declared extinct. Every biological split uses the same speciation process;
the artificial binary node representing a sampled ancestor is not a split.

Hidden lineage histories are marginalized numerically with a hidden-count
forward calculation. The likelihood does not substitute sampled-tree diversity
for total diversity. See the [likelihood derivation and implementation plan](doc/diversity-dependent-fbd/SPECIMEN_PLAN.md).

## Obtain and build this branch

The function is experimental and is not part of an unmodified RevBayes release.
With Git, a C++23 compiler, Meson, Ninja and Boost installed:

```sh
git clone --branch ddfbd-model https://github.com/sergeitarasov/revbayes.git
cd revbayes
projects/meson/generate.sh
meson setup build
ninja -C build
```

See the [build instructions](projects/meson/README.md) for dependencies and Boost
configuration. The development build used Boost 1.88.0. On Apple Clang/libc++,
if the bundled range-v3/meta headers conflict with standard-library forward
declarations, configure with
`meson setup build -Dcpp_args=-DMETA_NO_STD_FORWARD_DECLARATIONS`
(or `meson configure build -Dcpp_args=-DMETA_NO_STD_FORWARD_DECLARATIONS` for an
existing build). This compatibility flag was used in the tested macOS build.
Local binaries and downloaded build dependencies are not included in Git.

## Run morphology-only inference

From the repository root:

```sh
build/rb examples/DDFBDP/morphology.Rev
```

The [complete Rev script](examples/DDFBDP/morphology.Rev) loads the supplied bear
data, estimates the diversification/preservation parameters and morphology
model, and updates topology, node ages, occurrence ages and ancestor status.
Its essential model connection is:

```rev
tree ~ dnDDFBDP(originAge=origin, lambda0=lambda0, alpha=alpha,
                mu=mu, psi=psi, rho=1, taxa=taxa,
                condition="sampling", maxHiddenLineages=128)

morphology ~ dnPhyloCTMC(tree=tree, Q=fnJC(2), branchRates=clock,
                        type="Standard", coding="variable")
morphology.clamp(morpho)
moves.append(mvCollapseExpandFossilBranch(tree, origin, weight=6))
```

The complete script defines the parameters, their priors, all moves and monitors;
the excerpt above is not a standalone script. Outputs are written to
`output/DDFBDP/morphology.log` and `output/DDFBDP/morphology.trees`. The example
runs 1,000 generations for demonstration, not a convergence assessment.

The data comprise 22 labels (8 living, 14 fossil) and 62 binary characters,
with morphology missing for four labels. Input provenance and adaptations are
explained in the [example guide](examples/DDFBDP/README.md). Fossil age bounds
are uncertainty intervals for individual occurrences, not species ranges.
When adapting the example, change the origin prior and character model to suit
the data. `coding="variable"` assumes constant characters were excluded.

## Conditioning and numerical settings

| Argument | Meaning |
|---|---|
| `originAge` | Single-lineage origin, older than the tree root; may be estimated. |
| `rho` | Independent sampling probability for each lineage alive at present. |
| `condition="time"` | Origin-started density with no observation conditioning. |
| `condition="sampling"` | Condition on at least one fossil or living observation. |
| `condition="survival"` | Condition on at least one sampled living descendant. |
| `maxHiddenLineages=128` | Numerical hidden-count cutoff; increase it to check sensitivity. |
| `numericalTolerance=1e-12` | Propagation tolerance; does not bound cutoff error. |
| `initialTree` | Optional sampled tree matching the taxa and occurrence-age bounds. |

At `alpha=0`, time- and sampling-conditioned densities match `dnBDSTP(r=0)`.
There is a documented native survival-conditioning difference: native BDSTP
divides the sampled-extant survival probability by `rho`, whereas this prior
uses the probability itself. Thus at `alpha=0`,
`lnP_native - lnP_DDFBDP = ln(rho)` for survival conditioning. The tests check
that factor explicitly; with `rho=1`, the densities agree directly. The new
prior retains BDSTP's legacy log-factorial approximation for the fixed label
normalization. Built-in usage documentation is available with `?dnDDFBDP`.

## Validation and its limits

The [saved integration report](tests/test_DDFBDP/results/integration.json) and
[numerical report](tests/test_DDFBDP/results/numerical_reference.json) document:

- 26 standalone numerical checks against analytic constant-rate densities and
  an independent DD solver; errors around 2.4e-12 or smaller in the reported comparisons.
- 46 native BDSTP comparisons: 39 direct and 7 with the survival correction
  above; maximum error 2.1e-12.
- A 200-generation morphology-only diagnostic chain with cache checks enabled,
  including 212 accepted fossil collapse/expand proposals and ancestor counts
  from 0 to 10 among saved states.
- Hidden-count cutoff checks on five saved states; maximum log-density change
  from cutoff 128 to 256 was 1.3e-10.
- Constrained topology, checkpoint restoration and continued MCMC sampling.
- All 52 regression checks for the previous named-species kernel still passing.

Reproduce the tests with:

```sh
python3 tests/test_DDFBDP/numerical_reference.py
python3 tests/test_DDFBDP/run_integration.py build/rb \
    --output tests/test_DDFBDP/results/integration.json
```

The integration runner uses an isolated temporary directory and defaults to
200 generations. The standalone tests require Python 3 and a C++17 compiler;
set `CXX` if `clang++` is unavailable. The full application requires C++23.

These checks establish numerical agreement and working MCMC integration, not
parameter recovery, posterior convergence or simulation-based calibration.
At least two observation labels are required. Unconditional prior simulation,
root-age starts, count conditioning, removal and skyline rates are not implemented.
A generated MCMC starting tree is an initializer, not a draw from the DD prior.

## Relationship to the earlier implementation

`dnDDFBDP` is separate from `dnDDFBD`, the earlier named-species extended-tree
model. The new specimen prior needs neither persistent species identities nor
an observation-tree projection, and uses standard fossil collapse/expand moves.
Repeated occurrences can have distinct specimen labels, but are not constrained
to a shared named species. The earlier model and its recovery experiments remain
available through [DDFBD.md](DDFBD.md); those experiments do not validate recovery
under this new observation model. The existing `dnBDSTP` implementation is unchanged.
