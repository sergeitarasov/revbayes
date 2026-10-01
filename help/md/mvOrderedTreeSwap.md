## name
mvOrderedTreeSwap
## title
Swap ordered subtrees while preserving budding-species orientations
## description
Propose a topology or orientation update on a rooted time tree without changing node ages or node indices.
## details
Draw two distinct nonroot nodes uniformly from the fixed node set. Ancestor-descendant pairs and pairs incompatible with strictly positive branch durations are rejected without resampling. Otherwise, exchange the selected subtrees in their exact original child slots. Swapping siblings reverses their order, which changes the budding-species orientation when child 0 denotes ancestral continuation. Swapping nonsiblings changes topology while preserving all other child slots.

The proposal is self-inverse with log Hastings ratio zero. Undo restores exact child ordering and parent links. Species identity, fossil occurrences, origin and extinction constraints remain the responsibility of the tree distribution and observation-tree likelihood; incompatible proposals are rejected through their target probability.

This move does not update ages. Combine it with appropriate age moves. Unlike ordinary unordered-tree rearrangements, it treats child order as part of the stochastic state.
## authors
RevBayes development contributors
## see_also
fnSpeciesObservationTree, mvNodeTimeSlideUniform
## example
moves.append(mvOrderedTreeSwap(species_tree, weight=10.0))
## references
