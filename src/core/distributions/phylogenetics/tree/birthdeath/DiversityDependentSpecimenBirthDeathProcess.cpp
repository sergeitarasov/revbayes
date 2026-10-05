#include "DiversityDependentSpecimenBirthDeathProcess.h"
#include "StartingTreeSimulator.h"
#include "TreeUtilities.h"
#include "TopologyNode.h"
#include "Clade.h"
#include "RbException.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include <algorithm>
#include <cmath>
#include <set>
using namespace RevBayesCore;
DiversityDependentSpecimenBirthDeathProcess::DiversityDependentSpecimenBirthDeathProcess(
    const TypedDagNode<double>* o, const TypedDagNode<double>* l,
    const TypedDagNode<double>* a, const TypedDagNode<double>* m,
    const TypedDagNode<double>* p, const TypedDagNode<double>* r,
    const std::vector<Taxon>& t, Tree* initial, const std::string& c,
    const DDFBD::Settings& s, std::int64_t pr)
    : AbstractBirthDeathProcess(o,c,t,true,initial),lambda(l),alpha(a),mu(m),psi(p),rho(r),settings(s),precision(pr)
{
    if (t.size()<2) throw RbException("dnDDFBDP currently requires at least two observation labels");
    if (c!="time" && c!="sampling" && c!="survival") throw RbException("Invalid dnDDFBDP condition");
    if (s.maxHidden<1) throw RbException("dnDDFBDP requires at least one hidden-count state beyond zero");
    for (const DagNode* node : {l,a,m,p,r}) addParameter(node);
    if (initial) setValue(initial->clone());
    else redrawValue(SimulationCondition::MCMC);
}
DiversityDependentSpecimenBirthDeathProcess* DiversityDependentSpecimenBirthDeathProcess::clone() const
{ return new DiversityDependentSpecimenBirthDeathProcess(*this); }
void DiversityDependentSpecimenBirthDeathProcess::setValue(Tree* tree,bool force)
{
    Tree* initialized=TreeUtilities::startingTreeInitializer(*tree,taxa,precision);
    AbstractRootedTreeDistribution::setValue(initialized,force);
}
void DiversityDependentSpecimenBirthDeathProcess::redrawValue()
{ throw RbException("dnDDFBDP: unconditional prior simulation with fixed specimen data is not implemented"); }
void DiversityDependentSpecimenBirthDeathProcess::redrawValue(SimulationCondition c)
{
    if (c!=SimulationCondition::MCMC) { redrawValue(); return; }
    if (starting_tree) { setValue(starting_tree->clone()); return; }
    RbVector<Clade> constraints;
    Clade all(taxa); all.setAge(getOriginAge()-1e-5); constraints.push_back(all);
    StartingTreeSimulator simulator;
    Tree* tree=simulator.simulateTree(taxa,constraints);
    AbstractRootedTreeDistribution::setValue(tree);
}
double DiversityDependentSpecimenBirthDeathProcess::computeLnProbabilityDivergenceTimes() const
{ return computeLnProbabilityTimes(); } // Override base root-survival normalization.
double DiversityDependentSpecimenBirthDeathProcess::computeLnProbabilityTimes() const
{
    const double impossible=RbConstants::Double::neginf;
    const double origin=getOriginAge();
    DDFBD::Parameters p{lambda->getValue(),alpha->getValue(),mu->getValue(),psi->getValue(),rho->getValue()};
    for (double x : {origin,p.lambda0,p.alpha,p.mu,p.psi,p.rho})
        if (!std::isfinite(x) || x<0) return impossible;
    if (origin==0 || p.rho>1 || value->getNumberOfTips()!=taxa.size()) return impossible;
    std::set<std::string> expected,seen;
    for (const auto& taxon:taxa) expected.insert(taxon.getName());
    std::vector<DDFBD::SpecimenEvent> events;
    std::size_t extant=0;
    for (const TopologyNode* node:value->getNodes()) {
        const double age=node->getAge();
        if (!std::isfinite(age) || age<0 || age>origin) return impossible;
        if (node->isTip()) {
            if (!expected.count(node->getName()) || !seen.insert(node->getName()).second) return impossible;
            const auto& taxon=node->getTaxon();
            // Newick branch-length subtraction can round an exact occurrence age.
            // Allow floating-point noise, not age moves away from a fixed sample.
            if (taxon.getMinAge()==taxon.getMaxAge() &&
                std::fabs(age-taxon.getMinAge())>1e-10*std::max(1.0,age)) return impossible;
            if (age==0) {
                if (node->isSampledAncestorTip()) return impossible;
                ++extant;
            } else {
                if (!node->isSampledAncestorTip() && node->getBranchLength()<=0) return impossible;
                events.push_back({age,node->isSampledAncestorTip() ?
                    DDFBD::SpecimenEventType::AncestorSample : DDFBD::SpecimenEventType::TerminalSample});
            }
        } else {
            if (node->getNumberOfChildren()!=2 || age==0) return impossible;
            unsigned sa=0;
            for (const auto* child:node->getChildren()) {
                if (child->isSampledAncestorTip()) ++sa;
                else if (child->getAge()>=age) return impossible;
            }
            if (sa>1) return impossible;
            if (!sa) events.push_back({age,DDFBD::SpecimenEventType::Birth});
        }
    }
    if (condition=="survival" && extant==0) return impossible;
    std::stable_sort(events.begin(),events.end(),[](const auto& a,const auto& b){
        if (a.age!=b.age) return a.age>b.age;
        return static_cast<int>(a.type)<static_cast<int>(b.type);
    });
    try {
        double result=DDFBD::specimenLogLikelihood(origin,events,extant,p,settings);
        if (condition!="time" && std::isfinite(result)) {
            const double logObs=DDFBD::logObservationProbability(origin,settings.maxHidden+taxa.size(),condition=="sampling",p,settings);
            if (!std::isfinite(logObs)) return impossible;
            result-=logObs;
        }
        return result;
    } catch(const std::exception& error) { throw RbException(std::string("dnDDFBDP: ")+error.what()); }
}
double DiversityDependentSpecimenBirthDeathProcess::lnProbTreeShape() const
{
    const double n=value->getNumberOfTips(),a=value->getNumberOfSampledAncestors();
    // Match BDSTP's legacy factorial approximation as well as its convention.
    // This label-normalization constant is fixed throughout an MCMC analysis.
    return (n-a-1)*RbConstants::LN2-RbMath::lnFactorial(static_cast<int>(n));
}
double DiversityDependentSpecimenBirthDeathProcess::lnProbNumTaxa(size_t,double,double,bool) const
{ throw RbException("dnDDFBDP does not implement count conditioning"); }
double DiversityDependentSpecimenBirthDeathProcess::pSurvival(double start,double end) const
{
    const DDFBD::Parameters p{lambda->getValue(),alpha->getValue(),mu->getValue(),psi->getValue(),rho->getValue()};
    return std::exp(DDFBD::logObservationProbability(start-end,settings.maxHidden+taxa.size(),false,p,settings));
}
double DiversityDependentSpecimenBirthDeathProcess::simulateDivergenceTime(double,double) const
{ throw RbException("dnDDFBDP divergence simulation is not an independent-lineage operation"); }
std::vector<double> DiversityDependentSpecimenBirthDeathProcess::simulateDivergenceTimes(
    size_t n,double origin,double present,double minimum,bool alwaysReturn) const
{
    // dnConstrainedTopology requests arbitrary valid times to initialize MCMC.
    // Its validation/prior simulation path must not use this initialization law.
    if (!alwaysReturn) throw RbException("dnDDFBDP prior simulation is not implemented");
    const double lower=std::max(present,minimum);
    if (!std::isfinite(origin) || !(origin>lower))
        throw RbException("dnDDFBDP initializer requires an origin older than the samples");
    std::vector<double> times(n);
    for (double& t:times) t=lower+(0.01+0.98*GLOBAL_RNG->uniform01())*(origin-lower);
    std::sort(times.begin(),times.end());
    return times;
}
void DiversityDependentSpecimenBirthDeathProcess::swapParameterInternal(const DagNode* oldP,const DagNode* newP)
{
    if (oldP==lambda) lambda=static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP==alpha) alpha=static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP==mu) mu=static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP==psi) psi=static_cast<const TypedDagNode<double>*>(newP);
    else if (oldP==rho) rho=static_cast<const TypedDagNode<double>*>(newP);
    else AbstractBirthDeathProcess::swapParameterInternal(oldP,newP);
}
