#include "Dist_DDFBDP.h"
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
#include "RlTaxon.h"
#include "RevNullObject.h"

using namespace RevLanguage;

Dist_DDFBDP::Dist_DDFBDP() : TypedDistribution<TimeTree>() {}
Dist_DDFBDP* Dist_DDFBDP::clone(void) const
{ return new Dist_DDFBDP(*this); }
const std::string& Dist_DDFBDP::getClassType(void)
{ static std::string name = "Dist_DDFBDP"; return name; }
const TypeSpec& Dist_DDFBDP::getClassTypeSpec(void)
{
    static TypeSpec spec(getClassType(), new TypeSpec(TypedDistribution<TimeTree>::getClassTypeSpec()));
    return spec;
}
const TypeSpec& Dist_DDFBDP::getTypeSpec(void) const { return getClassTypeSpec(); }
std::string Dist_DDFBDP::getDistributionFunctionName(void) const
{ return "DiversityDependentSpecimenFBD"; }
std::vector<std::string> Dist_DDFBDP::getDistributionFunctionAliases(void) const
{ return {"DDFBDP"}; }

RevBayesCore::DiversityDependentSpecimenBirthDeathProcess*
Dist_DDFBDP::createDistribution(void) const
{
    RevBayesCore::DDFBD::Settings settings;
    settings.maxHidden = static_cast<const Natural&>(maxHidden->getRevObject()).getValue();
    settings.tolerance = static_cast<const RealPos&>(tolerance->getRevObject()).getValue();
    return new RevBayesCore::DiversityDependentSpecimenBirthDeathProcess(
        static_cast<const RealPos&>(origin->getRevObject()).getDagNode(),
        static_cast<const RealPos&>(lambda0->getRevObject()).getDagNode(),
        static_cast<const RealPos&>(alpha->getRevObject()).getDagNode(),
        static_cast<const RealPos&>(mu->getRevObject()).getDagNode(),
        static_cast<const RealPos&>(psi->getRevObject()).getDagNode(),
        static_cast<const Probability&>(rho->getRevObject()).getDagNode(),
        static_cast<const ModelVector<Taxon>&>(taxa->getRevObject()).getValue(),
        initialTree->getRevObject()==RevNullObject::getInstance() ? nullptr :
            static_cast<const TimeTree&>(initialTree->getRevObject()).getValue().clone(),
        static_cast<const RlString&>(condition->getRevObject()).getValue(), settings,
        static_cast<const Natural&>(agePrecision->getRevObject()).getValue());
}

const MemberRules& Dist_DDFBDP::getParameterRules(void) const
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
        rules.push_back(new ArgumentRule("taxa", ModelVector<Taxon>::getClassTypeSpec(), "Fossil specimens and living taxa, with occurrence-age bounds.", ArgumentRule::BY_VALUE, ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("initialTree", TimeTree::getClassTypeSpec(), "Optional sampled tree for initialization.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL));
        rules.push_back(new ArgumentRule("ageCheckPrecision", Natural::getClassTypeSpec(), "Precision used to check imported sample ages.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(4)));
        rules.push_back(new OptionRule("condition", new RlString("time"), {"time", "sampling", "survival"}, "Origin only, any observation, or sampled extant survival."));
        rules.push_back(new ArgumentRule("maxHiddenLineages", Natural::getClassTypeSpec(), "Killed hidden-count truncation; verify convergence by increasing it.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(128)));
        rules.push_back(new ArgumentRule("numericalTolerance", RealPos::getClassTypeSpec(), "Matrix exponential tolerance (not a hidden-count error bound).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1e-12)));
        initialized = true;
    }
    return rules;
}

void Dist_DDFBDP::setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var)
{
    if (name == "originAge") origin = var;
    else if (name == "lambda0") lambda0 = var;
    else if (name == "alpha") alpha = var;
    else if (name == "mu") mu = var;
    else if (name == "psi") psi = var;
    else if (name == "rho") rho = var;
    else if (name == "taxa") taxa = var;
    else if (name == "ageCheckPrecision") agePrecision = var;
    else if (name == "initialTree") initialTree = var;
    else if (name == "condition") condition = var;
    else if (name == "maxHiddenLineages") maxHidden = var;
    else if (name == "numericalTolerance") tolerance = var;
    else Distribution::setConstParameter(name,var);
}
