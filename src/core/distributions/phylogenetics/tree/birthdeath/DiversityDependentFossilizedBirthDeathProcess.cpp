#include "DiversityDependentFossilizedBirthDeathProcess.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

using namespace RevBayesCore;

DiversityDependentFossilizedBirthDeathProcess::DiversityDependentFossilizedBirthDeathProcess(
    const TypedDagNode<double>* origin, const TypedDagNode<double>* lambda0,
    const TypedDagNode<double>* alpha, const TypedDagNode<double>* mu,
    const TypedDagNode<double>* psi, const TypedDagNode<double>* rho,
    const std::vector<std::string>& occurrenceSpecies,
    const TypedDagNode<RbVector<double> >* occurrenceAges,
    const std::vector<std::string>& sampledExtant, const Tree& initialTree,
    const std::string& condition, const DDFBD::Settings& settings) :
    TypedDistribution<Tree>(initialTree.clone()), originAge(origin), birthRate(lambda0),
    feedback(alpha), deathRate(mu), fossilRate(psi), extantSampling(rho),
    fossilAges(occurrenceAges), fossilSpecies(occurrenceSpecies),
    sampledLiving(sampledExtant.begin(), sampledExtant.end()), initial(initialTree),
    conditioning(condition), numerical(settings)
{
    for (const DagNode* p : std::vector<const DagNode*>{origin,lambda0,alpha,mu,psi,rho,occurrenceAges})
        addParameter(p);
    if (conditioning != "none" && conditioning != "sampling" && conditioning != "sampledExtant")
        throw RbException("DDFBD condition must be none, sampling, or sampledExtant");
    if (sampledLiving.size() != sampledExtant.size())
        throw RbException("DDFBD duplicate sampledExtant species");
    if (fossilSpecies.size() != fossilAges->getValue().size())
        throw RbException("DDFBD occurrenceSpecies and occurrenceAges must have equal lengths");
    for (const TopologyNode* tip : initial.getNodes()) {
        if (!tip->isTip()) continue;
        if (!livingStatus.emplace(tip->getName(), tip->getAge() == 0).second)
            throw RbException("DDFBD initial tree contains duplicate species names");
    }
    for (const std::string& name : fossilSpecies)
        if (!livingStatus.count(name)) throw RbException("DDFBD unknown occurrence species: " + name);
    for (const std::string& name : sampledLiving)
        if (!livingStatus.count(name) || !livingStatus.at(name))
            throw RbException("DDFBD sampledExtant species must have a present-day endpoint: " + name);
    for (const auto& entry : livingStatus)
        if (!sampledLiving.count(entry.first) &&
            std::find(fossilSpecies.begin(), fossilSpecies.end(), entry.first) == fossilSpecies.end())
            throw RbException("DDFBD every retained species needs an occurrence or extant sample: " + entry.first);
    if (livingStatus.empty()) throw RbException("DDFBD requires a nonempty initialTree");
    if (conditioning == "sampledExtant" && sampledLiving.empty())
        throw RbException("DDFBD sampledExtant conditioning requires an extant sample");
}

DiversityDependentFossilizedBirthDeathProcess*
DiversityDependentFossilizedBirthDeathProcess::clone(void) const
{
    return new DiversityDependentFossilizedBirthDeathProcess(*this);
}

double DiversityDependentFossilizedBirthDeathProcess::computeLnProbability(void)
{
    const double impossible = RbConstants::Double::neginf;
    const double origin = originAge->getValue();
    if (!std::isfinite(origin) || origin <= 0 || value->getNumberOfTips() != livingStatus.size())
        return impossible;
    if (fossilAges->getValue().size() != fossilSpecies.size()) return impossible;
    std::map<std::string, double> attachment, endpoint, oldest, youngest;
    std::vector<DDFBD::Event> events;
    bool valid = true;
    std::function<std::string(const TopologyNode&)> visit = [&](const TopologyNode& node) {
        if (!std::isfinite(node.getAge()) || node.getAge() < 0 || node.getAge() >= origin || node.isSampledAncestorTip())
            valid = false;
        if (node.isTip()) {
            const std::string& name = node.getName();
            auto status = livingStatus.find(name);
            if (status == livingStatus.end() || endpoint.count(name)) valid = false;
            else if (status->second != (node.getAge() == 0)) valid = false;
            endpoint[name] = node.getAge();
            return name;
        }
        if (node.getNumberOfChildren() != 2) { valid = false; return std::string(); }
        for (std::size_t c = 0; c < 2; ++c)
            if (node.getChild(c).getAge() >= node.getAge()) valid = false;
        const std::string ancestor = visit(node.getChild(0));
        const std::string daughter = visit(node.getChild(1));
        attachment[daughter] = node.getAge();
        events.push_back({node.getAge(), DDFBD::EventType::Birth});
        return ancestor;
    };
    const std::string initialSpecies = visit(value->getRoot());
    attachment[initialSpecies] = origin;
    if (!valid || endpoint.size() != livingStatus.size()) return impossible;
    for (const auto& entry : endpoint) {
        oldest[entry.first] = 0;
        youngest[entry.first] = origin;
    }
    for (std::size_t i = 0; i < fossilSpecies.size(); ++i) {
        const double age = fossilAges->getValue()[i];
        const std::string& name = fossilSpecies[i];
        if (!std::isfinite(age) || age < 0 || age < endpoint.at(name) || age > attachment.at(name))
            return impossible;
        oldest[name] = std::max(oldest[name],age);
        youngest[name] = std::min(youngest[name],age);
    }
    std::size_t retainedUnsampledLiving = 0;
    for (const auto& entry : endpoint) {
        const std::string& name = entry.first;
        if (sampledLiving.count(name)) youngest[name] = 0;
        if (oldest[name] < entry.second || attachment.at(name) < oldest[name]) return impossible;
        events.push_back({oldest[name], DDFBD::EventType::Protect});
        if (entry.second > 0) events.push_back({entry.second, DDFBD::EventType::Death});
        else if (!sampledLiving.count(name)) ++retainedUnsampledLiving;
    }
    std::sort(events.begin(), events.end(), [](const DDFBD::Event& a, const DDFBD::Event& b) {
        if (a.age != b.age) return a.age > b.age;
        return static_cast<int>(a.type) < static_cast<int>(b.type);
    });
    const DDFBD::Parameters p{birthRate->getValue(), feedback->getValue(), deathRate->getValue(),
                              fossilRate->getValue(), extantSampling->getValue()};
    // Scalar proposals may temporarily leave their Rev type/prior support.
    // Those states have zero density, rather than being fatal configuration errors.
    for (double parameter : {p.lambda0,p.alpha,p.mu,p.psi,p.rho})
        if (!std::isfinite(parameter) || parameter < 0) return impossible;
    if (p.rho > 1) return impossible;
    try {
        double logDensity = DDFBD::logLikelihood(origin, events, fossilSpecies.size(),
            sampledLiving.size(), retainedUnsampledLiving, p, numerical);
        // Orientations remain explicit in the state; do NOT multiply by 2^(n-1).
        logDensity -= std::lgamma(static_cast<double>(livingStatus.size())+1);
        if (conditioning != "none" && std::isfinite(logDensity)) {
            const double logSample = DDFBD::logObservationProbability(origin,
                numerical.maxHidden+livingStatus.size(), conditioning == "sampling", p, numerical);
            if (logSample == impossible) return impossible;
            logDensity -= logSample;
        }
        return logDensity;
    } catch (const std::invalid_argument& error) {
        // Invalid schedule here is an implementation/data-contract error, not a
        // silently rejected MCMC proposal. Tree support was checked above.
        throw RbException(std::string("DDFBD: ")+error.what());
    } catch (const std::runtime_error& error) {
        throw RbException(std::string("DDFBD numerical failure: ")+error.what());
    }
}

void DiversityDependentFossilizedBirthDeathProcess::redrawValue(void)
{
    throw RbException("DDFBD prior redraw with fixed named occurrence assignments is not implemented; use the complete-history simulator for validation");
}

void DiversityDependentFossilizedBirthDeathProcess::redrawValue(SimulationCondition condition)
{
    if (condition == SimulationCondition::MCMC) {
        // Valid starting state only, not a claimed draw from the prior.
        setValue(initial.clone());
    } else redrawValue();
}

void DiversityDependentFossilizedBirthDeathProcess::swapParameterInternal(const DagNode* oldP, const DagNode* newP)
{
    if (oldP == originAge) originAge = static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP == birthRate) birthRate = static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP == feedback) feedback = static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP == deathRate) deathRate = static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP == fossilRate) fossilRate = static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP == extantSampling) extantSampling = static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP == fossilAges) fossilAges = static_cast<const TypedDagNode<RbVector<double> >*>(newP);
}
