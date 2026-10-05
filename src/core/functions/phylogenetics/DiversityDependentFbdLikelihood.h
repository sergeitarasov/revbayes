#ifndef RB_DIVERSITY_DEPENDENT_FBD_LIKELIHOOD_H
#define RB_DIVERSITY_DEPENDENT_FBD_LIKELIHOOD_H

#include <cstddef>
#include <vector>

namespace RevBayesCore { namespace DDFBD {

/** Numerical kernels for DD tree densities. No DAG dependencies.
 * Hidden counts are truncated with a killed, never reflecting, boundary.
 * See doc/diversity-dependent-fbd/EXTENDED_TREE_DERIVATION.md for the measure.
 */
enum class BirthModel { Exponential, DDDLinear, DDDPower };
struct Parameters {
    double lambda0;
    double alpha;
    double mu;
    double psi;
    double rho;
    BirthModel birthModel = BirthModel::Exponential;
    double capacity = 1;
    double birth(std::size_t total) const;
};

enum class EventType { Birth, Protect, Death };
struct Event {
    double age;
    EventType type;
};

struct Settings {
    std::size_t maxHidden = 128;
    double tolerance = 1e-12;
};

// Specimen observations: an ancestral sample continues the retained skeleton;
// a terminal sample demotes its still-living continuation into the hidden count.
enum class SpecimenEventType { Birth, AncestorSample, TerminalSample };
struct SpecimenEvent { double age; SpecimenEventType type; };
// Raw density, before the labelled non-oriented sampled-tree shape factor and
// observation conditioning. Events must be oldest first; times are positive.
// The origin starts with one retained lineage and no hidden lineages.
double specimenLogLikelihood(double origin, const std::vector<SpecimenEvent>& events,
    std::size_t extantSamples, const Parameters& parameters, const Settings& settings);

/** Raw oriented density; excludes label factorial and sampling conditioning.
 * events must be oldest first, ties Birth then Protect then Death. Every species
 * has one Protect event (its oldest observation) and either a Death endpoint or
 * is counted in sampledLiving/unsampledLiving. Fossil records contribute psi^M.
 * Throws on malformed schedules/numerical failure; returns -infinity for zero
 * density due to valid parameter boundaries. Does not change numerical settings.
 */
double logLikelihood(double origin, const std::vector<Event>& events,
                     std::size_t fossilRecords, std::size_t sampledLiving,
                     std::size_t unsampledLiving, const Parameters& parameters,
                     const Settings& settings);

/** Complete extant tree likelihood, in DDD's phylogeny (1) or branching-time
 * (0) convention. Crown initialization excludes the root birth factor; stem
 * initialization includes every branching age. Conditions: 0 age only,
 * 1 survival, 2 extant count. No fossils, missing extant species or cond=3.
 * branchingAges contains all internal ages, oldest first (including crown).
 */
double extantTreeLogLikelihood(double startAge, const std::vector<double>& branchingAges,
    bool crown, unsigned condition, bool branchingTimes, const Parameters& parameters,
    const Settings& settings);

/** Probability of no observations from one initial species, evaluated in a
 * count process (not the embedding kernel), with a killed total-count boundary.
 * fossilSampling=true includes fossil observations in the conditioning event.
 */
double noObservationProbability(double origin, std::size_t maxTotal,
                                bool fossilSampling, const Parameters& parameters,
                                const Settings& settings);

/** Direct log P(at least one observation), with an absorbing observed state.
 * Unlike 1-noObservationProbability, neither cutoff leakage nor cancellation
 * is mistaken for sampling. Use this for ascertainment normalization.
 */
double logObservationProbability(double origin, std::size_t maxTotal,
                                 bool fossilSampling, const Parameters& parameters,
                                 const Settings& settings);

} }
#endif
