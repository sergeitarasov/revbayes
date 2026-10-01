## name
fnSpeciesObservationTree
## title
Project an oriented budding-species tree onto dated character samples
## description
Construct a deterministic TimeTree for molecular or morphological data whose observations occur at specified times along named species lineages.
## details
The input tree has one terminal endpoint per named species. At every bifurcation, child 0 continues the ancestral species and child 1 starts the new species. Endpoint names identify species; endpoint ages are extinction times or zero for living species. This orientation is part of the model and must be preserved by compatible tree proposals.

The equally sized vectors sampleNames, speciesNames and ages identify character samples, their species, and their actual ages before the present. Sample names must be unique. Repeated species names are allowed. Fossil observations with surviving sampled descendants become zero-length sampled-ancestor tips; an observation with no sampled descendants becomes a terminal tip at its observation age. Unsampled species endpoints and resulting unary branching nodes are pruned. Tip indices follow sampleNames order.

The initial species can have samples on its stem older than the species-tree root. The upstream species-tree distribution must constrain these to be younger than the origin. Other samples must occur after their species birth and before their endpoint. Out-of-support dynamic ages or placements raise MATH_ERROR, allowing Metropolis-Hastings rejection; malformed static identifiers raise an ordinary error.

This function supplies the character-data genealogy only. It does not add a fossil-sampling likelihood. Fossil occurrences and living samples must also be represented in the species-tree distribution's sampling data. A species-level composite character row must have an explicitly chosen observation time; do not duplicate it for every occurrence.

## authors
RevBayes development contributors
## see_also
fnPruneTree, dnPhyloCTMC
## example
# species_tree is an oriented extended species tree with endpoints A and B.
# A has a historical character sample and a living sample; B has one fossil sample.
observations := fnSpeciesObservationTree(species_tree,
    sampleNames=v("A_fossil", "A_living", "B_fossil"),
    speciesNames=v("A", "A", "B"), ages=v(3.0, 0.0, 2.0))
# Use observations as tree= in dnPhyloCTMC; data row names must match sampleNames.
## references
Stadler T, Gavryushkina A, Warnock RCM, Drummond AJ, Heath TA. 2018. The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Journal of Theoretical Biology.
