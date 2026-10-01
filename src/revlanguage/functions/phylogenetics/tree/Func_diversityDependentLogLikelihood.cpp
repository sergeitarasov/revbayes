#include "Func_diversityDependentLogLikelihood.h"
#include "DiversityDependentTreeLikelihoodFunction.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlTimeTree.h"
#include "RlString.h"
#include "RealPos.h"
#include "Natural.h"
#include "OptionRule.h"
using namespace RevLanguage;
Func_diversityDependentLogLikelihood::Func_diversityDependentLogLikelihood() : TypedFunction<Real>() {}
Func_diversityDependentLogLikelihood* Func_diversityDependentLogLikelihood::clone() const
{ return new Func_diversityDependentLogLikelihood(*this); }
const std::string& Func_diversityDependentLogLikelihood::getClassType()
{ static std::string s="Func_diversityDependentLogLikelihood"; return s; }
const TypeSpec& Func_diversityDependentLogLikelihood::getClassTypeSpec()
{ static TypeSpec s(getClassType(),new TypeSpec(Function::getClassTypeSpec())); return s; }
const TypeSpec& Func_diversityDependentLogLikelihood::getTypeSpec() const { return getClassTypeSpec(); }
std::string Func_diversityDependentLogLikelihood::getFunctionName() const
{ return "fnDiversityDependentLogLikelihood"; }
RevBayesCore::TypedFunction<double>* Func_diversityDependentLogLikelihood::createFunction() const
{
    auto number = [this](size_t i) { return static_cast<const RealPos&>(args[i].getVariable()->getRevObject()).getDagNode(); };
    auto text = [this](size_t i) { return static_cast<const RlString&>(args[i].getVariable()->getRevObject()).getValue(); };
    auto model = RevBayesCore::DDFBD::BirthModel::Exponential;
    if (text(6)=="DDDlinear") model=RevBayesCore::DDFBD::BirthModel::DDDLinear;
    if (text(6)=="DDDpower") model=RevBayesCore::DDFBD::BirthModel::DDDPower;
    unsigned condition = text(8)=="none" ? 0 : (text(8)=="survival" ? 1 : 2);
    RevBayesCore::DDFBD::Settings settings;
    settings.maxHidden=static_cast<const Natural&>(args[10].getVariable()->getRevObject()).getValue();
    settings.tolerance=static_cast<const RealPos&>(args[11].getVariable()->getRevObject()).getValue();
    return new RevBayesCore::DiversityDependentTreeLikelihoodFunction(
        static_cast<const TimeTree&>(args[0].getVariable()->getRevObject()).getDagNode(),
        number(1),number(2),number(3),number(4),number(5),model,text(7)=="crown",
        condition,text(9)=="branchingTimes",settings);
}
const ArgumentRules& Func_diversityDependentLogLikelihood::getArgumentRules() const
{
    static ArgumentRules rules;
    static bool initialized=false;
    if (!initialized) {
        rules.push_back(new ArgumentRule("tree",TimeTree::getClassTypeSpec(),"Complete extant time tree.",ArgumentRule::BY_CONSTANT_REFERENCE,ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("lambda0",RealPos::getClassTypeSpec(),"Speciation scale.",ArgumentRule::BY_CONSTANT_REFERENCE,ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("mu",RealPos::getClassTypeSpec(),"Constant extinction rate.",ArgumentRule::BY_CONSTANT_REFERENCE,ArgumentRule::ANY));
        rules.push_back(new ArgumentRule("alpha",RealPos::getClassTypeSpec(),"Exponential decline coefficient.",ArgumentRule::BY_CONSTANT_REFERENCE,ArgumentRule::ANY,new RealPos(0.0)));
        rules.push_back(new ArgumentRule("K",RealPos::getClassTypeSpec(),"DDD equilibrium diversity.",ArgumentRule::BY_CONSTANT_REFERENCE,ArgumentRule::ANY,new RealPos(10.0)));
        rules.push_back(new ArgumentRule("originAge",RealPos::getClassTypeSpec(),"Required for stem start; ignored for crown.",ArgumentRule::BY_CONSTANT_REFERENCE,ArgumentRule::ANY,new RealPos(0.0)));
        rules.push_back(new OptionRule("rateModel",new RlString("exponential"),{"exponential","DDDlinear","DDDpower"},"Speciation rate law."));
        rules.push_back(new OptionRule("start",new RlString("crown"),{"crown","stem"},"Initial lineage convention."));
        rules.push_back(new OptionRule("condition",new RlString("survival"),{"none","survival","nTaxa"},"DDD ascertainment convention."));
        rules.push_back(new OptionRule("density",new RlString("DDDphylogeny"),{"DDDphylogeny","branchingTimes"},"Likelihood measure."));
        rules.push_back(new ArgumentRule("maxHiddenLineages",Natural::getClassTypeSpec(),"Killed hidden-count cutoff.",ArgumentRule::BY_VALUE,ArgumentRule::CONSTANT,new Natural(128)));
        rules.push_back(new ArgumentRule("numericalTolerance",RealPos::getClassTypeSpec(),"Uniformization tolerance.",ArgumentRule::BY_VALUE,ArgumentRule::CONSTANT,new RealPos(1e-12)));
        initialized=true;
    }
    return rules;
}
