#include "Dist_diversityDependentFBD.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlString.h"
#include "RevVariable.h"
#include "TypeSpec.h"

using namespace RevLanguage;

Dist_diversityDependentFBD::Dist_diversityDependentFBD() : TypedDistribution<TimeTree>() {}
Dist_diversityDependentFBD* Dist_diversityDependentFBD::clone(void) const
{ return new Dist_diversityDependentFBD(*this); }
const std::string& Dist_diversityDependentFBD::getClassType(void)
{ static std::string name = "Dist_diversityDependentFBD"; return name; }
const TypeSpec& Dist_diversityDependentFBD::getClassTypeSpec(void)
{
    static TypeSpec spec(getClassType(), new TypeSpec(TypedDistribution<TimeTree>::getClassTypeSpec()));
    return spec;
}
const TypeSpec& Dist_diversityDependentFBD::getTypeSpec(void) const { return getClassTypeSpec(); }
std::string Dist_diversityDependentFBD::getDistributionFunctionName(void) const
{ return "DiversityDependentFBD"; }
std::vector<std::string> Dist_diversityDependentFBD::getDistributionFunctionAliases(void) const
{ return {"DDFBD"}; }

RevBayesCore::DiversityDependentFossilizedBirthDeathProcess*
Dist_diversityDependentFBD::createDistribution(void) const
{
    RevBayesCore::DDFBD::Settings settings;
    settings.maxHidden = static_cast<const Natural&>(maxHidden->getRevObject()).getValue();
    settings.tolerance = static_cast<const RealPos&>(tolerance->getRevObject()).getValue();
    return new RevBayesCore::DiversityDependentFossilizedBirthDeathProcess(
        static_cast<const RealPos&>(origin->getRevObject()).getDagNode(),
        static_cast<const RealPos&>(lambda0->getRevObject()).getDagNode(),
        static_cast<const RealPos&>(alpha->getRevObject()).getDagNode(),
        static_cast<const RealPos&>(mu->getRevObject()).getDagNode(),
        static_cast<const RealPos&>(psi->getRevObject()).getDagNode(),
        static_cast<const Probability&>(rho->getRevObject()).getDagNode(),
        static_cast<const ModelVector<RlString>&>(occurrenceSpecies->getRevObject()).getValue(),
        static_cast<const ModelVector<RealPos>&>(occurrenceAges->getRevObject()).getDagNode(),
        static_cast<const ModelVector<RlString>&>(sampledExtant->getRevObject()).getValue(),
        static_cast<const TimeTree&>(initialTree->getRevObject()).getValue(),
        static_cast<const RlString&>(condition->getRevObject()).getValue(), settings);
}

const MemberRules& Dist_diversityDependentFBD::getParameterRules(void) const
{
    static MemberRules rules;
    static bool initialized = false;
    if (!initialized) {
        rules.push_back(new ArgumentRule("originAge", RealPos::getClassTypeSpec(), "Age of the single-species origin.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("lambda0", RealPos::getClassTypeSpec(), "Per-species birth rate at diversity one.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("alpha", RealPos::getClassTypeSpec(), "Nonnegative diversity feedback in exp(-alpha*(N-1)).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("mu", RealPos::getClassTypeSpec(), "Constant per-species extinction rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("psi", RealPos::getClassTypeSpec(), "Constant per-species fossil occurrence rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("rho", Probability::getClassTypeSpec(), "Independent extant sampling probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("occurrenceSpecies", ModelVector<RlString>::getClassTypeSpec(), "Assigned species of every supplied fossil record, including repeats.", ArgumentRule::BY_VALUE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("occurrenceAges", ModelVector<RealPos>::getClassTypeSpec(), "True ages of individual fossil records; these may be stochastic.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("sampledExtant", ModelVector<RlString>::getClassTypeSpec(), "Names of species sampled alive at present.", ArgumentRule::BY_VALUE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("initialTree", TimeTree::getClassTypeSpec(), "Oriented extended tree: child zero continues the ancestor, tips are extinction/present endpoints. Fixes living/extinct status.", ArgumentRule::BY_VALUE, ArgumentRule::ANY));
        rules.push_back(new OptionRule("condition", new RlString("none"), {"none", "sampling", "sampledExtant"}, "Ascertainment event; origin is not survival conditioning."));
        rules.push_back(new ArgumentRule("maxHiddenLineages", Natural::getClassTypeSpec(), "Killed hidden-count truncation; verify convergence by increasing it.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(128)));
        rules.push_back(new ArgumentRule("numericalTolerance", RealPos::getClassTypeSpec(), "Matrix exponential tolerance (not a hidden-count error bound).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1e-12)));
        initialized = true;
    }
    return rules;
}

void Dist_diversityDependentFBD::setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var)
{
    if (name == "originAge") origin = var;
    else if (name == "lambda0") lambda0 = var;
    else if (name == "alpha") alpha = var;
    else if (name == "mu") mu = var;
    else if (name == "psi") psi = var;
    else if (name == "rho") rho = var;
    else if (name == "occurrenceSpecies") occurrenceSpecies = var;
    else if (name == "occurrenceAges") occurrenceAges = var;
    else if (name == "sampledExtant") sampledExtant = var;
    else if (name == "initialTree") initialTree = var;
    else if (name == "condition") condition = var;
    else if (name == "maxHiddenLineages") maxHidden = var;
    else if (name == "numericalTolerance") tolerance = var;
    else Distribution::setConstParameter(name,var);
}
