# Independent DD-FBD simulation checks

`simulate_history.py` generates complete budding histories using Gillespie total
population hazards. It retains species identities, records every observed fossil,
and samples living species independently at present. The independent history
evaluator uses identity-specific event densities. The pruner retains a terminal
endpoint at the actual extinction/present time of every observed named species;
the first child at each retained split is the ancestral continuation A.

Run from the repository root:

```sh
python3 validation/DDFBD/simulate_history.py --self-test
python3 validation/DDFBD/simulate_history.py --seed 7 --output /tmp/ddfbd-seed7
python3 tests/test_DDFBD/numerical_reference.py
```

The committed `example_seed7` is one unconditional simulation, not a selected
recovery example. It has six observed species and nine fossils. Files are:

- `history.json`: all simulated species events, including hidden lineages;
- `extended_tree.json`: exact node ages, retained named species, and kernel events;
- `extended_tree.nwk`: oriented extended tree, with A first and D second;
- `extended_tree_annotated.nwk`: the same tree with explicit node-age annotations;
- `occurrences.tsv`: complete observed fossil records, with stable specimen IDs;
- `species.tsv`: endpoint ages, sampling indicators, and both structural attachment
  and true biological origination ages;
- `manifest.json`: parameters, seed, and generation conventions.

The Newick root branch denotes the origin stem. Generic Newick readers may ignore
that branch and may shift fossil-only trees because they lack a present-day tip;
use the JSON node ages and the species endpoint table when constructing a Rev
TimeTree. Never reinterpret structural attachment as biological species birth:
unsampled ancestors can be suppressed before a species' oldest occurrence.

An unconditional draw can have no observed species. In that case the JSON and
empty-header tables are written, and no Newick is produced; the simulator does
not silently redraw. Use a new output directory for each draw.

These checks validate event-density calculations and pruning invariants. They
are not simulation-based calibration or evidence of parameter recovery.
