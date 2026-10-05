#include "DiversityDependentFbdLikelihood.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace {
const double negativeInfinity = -std::numeric_limits<double>::infinity();

void validate(const RevBayesCore::DDFBD::Parameters& p,
              const RevBayesCore::DDFBD::Settings& settings, double origin)
{
    for (double x : {p.lambda0, p.alpha, p.mu, p.psi, p.rho, origin})
        if (!std::isfinite(x) || x < 0)
            throw std::invalid_argument("DDFBD requires finite nonnegative rates and ages");
    if (p.birthModel != RevBayesCore::DDFBD::BirthModel::Exponential &&
        (!std::isfinite(p.capacity) || p.capacity <= 0 || p.lambda0 <= p.mu || p.mu == 0))
        throw std::invalid_argument("DDD rate models require lambda0 > mu > 0 and finite K > 0");
    if (p.rho > 1 || origin == 0)
        throw std::invalid_argument("DDFBD requires rho in [0,1] and a positive origin");
    if (!std::isfinite(settings.tolerance) || settings.tolerance < 1e-15 || settings.tolerance > 1e-4)
        throw std::invalid_argument("DDFBD tolerance must be between 1e-15 and 1e-4");
    if (settings.maxHidden > 1000000)
        throw std::invalid_argument("DDFBD hidden-count limit exceeds supported allocation");
}

double normalize(std::vector<double>& v)
{
    const double scale = *std::max_element(v.begin(), v.end());
    if (!std::isfinite(scale))
        throw std::runtime_error("DDFBD numerical overflow");
    if (scale == 0) return negativeInfinity;
    for (double& x : v) {
        if (!std::isfinite(x) || x < 0)
            throw std::runtime_error("DDFBD lost positivity");
        x /= scale;
    }
    return std::log(scale);
}

// Positive uniformization of a Metzler matrix. This is NOT assumed stochastic:
// B=I+A/nu can have row/column sums greater than one. The infinity norm bounds
// the Taylor-series remainder. Substeps keep nu*dt <= 4 and prevent overflow.
double propagate(std::vector<double>& v, double duration,
                 const std::vector<double>& diagonal,
                 const std::vector<double>& up, const std::vector<double>& down,
                 double tolerance, const std::vector<double>& toObserved = {})
{
    if (duration == 0) return 0;
    const std::size_t size = v.size();
    double nu = 0;
    for (double d : diagonal) {
        if (!std::isfinite(d) || d > 0)
            throw std::runtime_error("DDFBD invalid diagonal hazard");
        nu = std::max(nu, -d);
    }
    if (nu == 0) return 0;
    const double work = std::ceil(nu * duration / 4.0);
    if (!std::isfinite(work) || work > 1000000)
        throw std::runtime_error("DDFBD propagation exceeds work limit; check rates, ages and cap");
    const std::size_t steps = std::max<std::size_t>(1, static_cast<std::size_t>(work));
    const double x = nu * (duration / steps);
    std::vector<double> bdiag(size);
    double bnorm = 0;
    for (std::size_t h = 0; h < size; ++h) {
        bdiag[h] = std::max(0.0, 1 + diagonal[h]/nu);
        if (toObserved.empty()) {
            const double rowsum = bdiag[h] + (h ? up[h-1]/nu : 0)
                                  + (h+1 < size ? down[h+1]/nu : 0);
            bnorm = std::max(bnorm, rowsum);
        } else {
            // The augmented conditioning process is substochastic. Its 1-norm
            // is well behaved even when many count states feed the absorber.
            const double colsum = bdiag[h] + (h ? down[h]/nu : 0)
                + (h+1 < size ? up[h]/nu : 0) + toObserved[h]/nu;
            bnorm = std::max(bnorm,colsum);
        }
    }
    double logScale = 0;
    std::vector<double> term(size), next(size), sum(size);
    for (std::size_t step = 0; step < steps; ++step) {
        term = v;
        sum = v;
        bool converged = false;
        for (std::size_t j = 1; j < 10000; ++j) {
            for (std::size_t h = 0; h < size; ++h) {
                next[h] = bdiag[h]*term[h];
                if (h) next[h] += up[h-1]*term[h-1]/nu;
                if (h+1 < size) next[h] += down[h+1]*term[h+1]/nu;
                next[h] *= x/j;
                sum[h] += next[h];
            }
            if (!toObserved.empty()) {
                double flux = 0;
                for (std::size_t h = 0; h+1 < size; ++h)
                    flux += toObserved[h]*term[h]/nu;
                flux *= x/j;
                next.back() += flux;
                sum.back() += flux;
            }
            term.swap(next);
            const double ratio = x * bnorm / (j+1);
            if (ratio < 1) {
                const double termNorm = toObserved.empty() ? *std::max_element(term.begin(), term.end())
                    : std::accumulate(term.begin(),term.end(),0.0);
                const double tail = termNorm * ratio / (1-ratio);
                const double scale = toObserved.empty() ? *std::max_element(sum.begin(), sum.end())
                    : std::accumulate(sum.begin(),sum.end(),0.0);
                if (tail <= (tolerance/steps)*scale) {
                    converged = true;
                    break;
                }
            }
        }
        if (!converged)
            throw std::runtime_error("DDFBD matrix exponential failed to converge");
        v.swap(sum);
        const double change = normalize(v);
        if (change == negativeInfinity) return change;
        logScale += change-x;
    }
    return logScale;
}

double interval(std::vector<double>& v, double duration, std::size_t k,
                std::size_t protectedCount, const RevBayesCore::DDFBD::Parameters& p,
                const RevBayesCore::DDFBD::Settings& settings)
{
    std::vector<double> diagonal(v.size()), up(v.size()), down(v.size());
    for (std::size_t h = 0; h < v.size(); ++h) {
        const double lambda = p.birth(k+h);
        diagonal[h] = -static_cast<double>(k+h)*(lambda+p.mu+p.psi);
        up[h] = (h+2*k-protectedCount)*lambda;
        down[h] = h*p.mu;
    }
    return propagate(v, duration, diagonal, up, down, settings.tolerance);
}

double logPower(double probability, std::size_t count)
{
    if (!count) return 0;
    return probability == 0 ? negativeInfinity : count*std::log(probability);
}
}

double RevBayesCore::DDFBD::Parameters::birth(std::size_t n) const
{
    if (!n) return 0;
    if (birthModel == BirthModel::DDDLinear)
        return std::max(0.0, lambda0-(lambda0-mu)*n/capacity);
    if (birthModel == BirthModel::DDDPower)
        return lambda0*std::pow(n+1.0, -std::log(lambda0/mu)/std::log1p(capacity));
    return lambda0*std::exp(-alpha*(n-1));
}

double RevBayesCore::DDFBD::logLikelihood(double origin, const std::vector<Event>& events,
    std::size_t fossilRecords, std::size_t sampledLiving, std::size_t unsampledLiving,
    const Parameters& p, const Settings& settings)
{
    validate(p, settings, origin);
    // With no extinction and complete extant sampling, every hidden lineage
    // would still be alive and sampled. Such histories have exactly zero weight.
    const std::size_t effectiveCap = p.mu == 0 && p.rho == 1 ? 0 : settings.maxHidden;
    std::vector<double> v(effectiveCap+1, 0);
    v[0] = 1;
    std::size_t k = 1, s = 0;
    double age = origin;
    double result = logPower(p.psi, fossilRecords);
    if (result == negativeInfinity) return result;
    for (const Event& event : events) {
        if (!std::isfinite(event.age) || event.age < 0 || event.age > age)
            throw std::invalid_argument("DDFBD events must be ordered oldest first within the origin");
        result += interval(v, age-event.age, k, s, p, settings);
        if (event.type == EventType::Birth) {
            if (!k) throw std::invalid_argument("DDFBD birth after all retained species ended");
            if (p.lambda0 == 0) return negativeInfinity;
            // Factor out the common rate in log space; exp(-alpha*(k-1))
            // can underflow even though its log density is finite.
            const double base = p.birth(k);
            if (p.birthModel == BirthModel::Exponential) {
                result += std::log(p.lambda0)-p.alpha*(k-1);
                for (std::size_t h = 0; h < v.size(); ++h) v[h] *= std::exp(-p.alpha*h);
            } else {
                if (base == 0) return negativeInfinity;
                result += std::log(base);
                for (std::size_t h = 0; h < v.size(); ++h) v[h] *= p.birth(k+h)/base;
            }
            ++k;
        } else if (event.type == EventType::Protect) {
            if (++s > k) throw std::invalid_argument("DDFBD protected count exceeds retained diversity");
        } else {
            if (!k || !s || event.age == 0)
                throw std::invalid_argument("DDFBD invalid extinction endpoint");
            result += logPower(p.mu, 1);
            --k;
            --s;
        }
        const double change = normalize(v);
        if (change == negativeInfinity || result == negativeInfinity) return negativeInfinity;
        result += change;
        age = event.age;
    }
    if (k != sampledLiving+unsampledLiving || s != k)
        throw std::invalid_argument("DDFBD endpoint/protected counts do not match species records");
    result += interval(v, age, k, s, p, settings);
    double endpoint = 0;
    double missingWeight = 1;
    for (double x : v) {
        endpoint += x*missingWeight;
        missingWeight *= 1-p.rho;
    }
    if (endpoint == 0)
        throw std::runtime_error("DDFBD endpoint density underflow; shorten time units or use a log-state solver");
    return result + std::log(endpoint) + logPower(p.rho, sampledLiving)
           + logPower(1-p.rho, unsampledLiving);
}

double RevBayesCore::DDFBD::noObservationProbability(double origin, std::size_t maxTotal,
    bool fossilSampling, const Parameters& p, const Settings& settings)
{
    validate(p, settings, origin);
    if (!maxTotal || maxTotal > 1000000)
        throw std::invalid_argument("DDFBD invalid conditioning count limit");
    std::vector<double> v(maxTotal+1,0), diag(v.size()), up(v.size()), down(v.size());
    v[1] = 1;
    for (std::size_t n = 0; n < v.size(); ++n) {
        up[n] = n*p.birth(n);
        down[n] = n*p.mu;
        diag[n] = -up[n]-down[n]-(fossilSampling ? n*p.psi : 0);
    }
    const double logScale = propagate(v, origin, diag, up, down, settings.tolerance);
    double weight = 1, result = 0;
    for (double x : v) {
        result += x*weight;
        weight *= 1-p.rho;
    }
    result *= std::exp(logScale);
    if (!std::isfinite(result) || result < 0 || result > 1+1e-10)
        throw std::runtime_error("DDFBD invalid conditioning probability");
    return std::min(1.0, result);
}

double RevBayesCore::DDFBD::logObservationProbability(double origin, std::size_t maxTotal,
    bool fossilSampling, const Parameters& p, const Settings& settings)
{
    validate(p,settings,origin);
    if (!maxTotal || maxTotal > 1000000)
        throw std::invalid_argument("DDFBD invalid conditioning count limit");
    if (p.rho == 0 && (!fossilSampling || p.psi == 0)) return negativeInfinity;
    if (p.lambda0 == 0) {
        // One initial species: a fossil before death OR an extant sample with
        // no earlier fossil. Evaluate the disjoint terms entirely in log space.
        const double psi = fossilSampling ? p.psi : 0;
        const double hazard = p.mu+psi;
        if (!std::isfinite(hazard)) throw std::runtime_error("DDFBD total hazard overflow");
        const double exposure = hazard*origin;
        const double logEvent = exposure == 0 ? std::log(hazard)+std::log(origin)
            : std::log(-std::expm1(-exposure));
        const double fossil = psi == 0 ? negativeInfinity :
            std::log(psi)-std::log(hazard)+logEvent;
        const double extant = p.rho == 0 ? negativeInfinity : std::log(p.rho)-hazard*origin;
        const double largest = std::max(fossil,extant);
        if (largest == negativeInfinity)
            throw std::runtime_error("DDFBD log observation probability overflow");
        return std::min(0.0, largest+std::log1p(std::exp(std::min(fossil,extant)-largest)));
    }
    // Last state records that an observation has occurred, irrespective of
    // subsequent genealogy. Cutoff exits are lost and never enter this state.
    std::vector<double> v(maxTotal+2,0), diag(v.size()), up(v.size()), down(v.size()), observed(v.size());
    v[1] = 1;
    for (std::size_t n = 0; n <= maxTotal; ++n) {
        const double births = n*p.birth(n);
        up[n] = n < maxTotal ? births : 0;
        down[n] = n*p.mu;
        observed[n] = fossilSampling ? n*p.psi : 0;
        diag[n] = -births-down[n]-observed[n];
    }
    const double logScale = propagate(v,origin,diag,up,down,settings.tolerance,observed);
    double result = v.back();
    for (std::size_t n = 1; n <= maxTotal; ++n) {
        const double sampleWeight = p.rho == 1 ? 1 : -std::expm1(n*std::log1p(-p.rho));
        result += v[n]*sampleWeight;
    }
    if (result == 0)
        throw std::runtime_error("DDFBD positive observation probability underflow; use a log-state solver");
    const double logResult = std::log(result)+logScale;
    if (!std::isfinite(logResult) || logResult > 1e-10)
        throw std::runtime_error("DDFBD invalid observation probability");
    return std::min(0.0,logResult);
}


double RevBayesCore::DDFBD::extantTreeLogLikelihood(double startAge,
    const std::vector<double>& branchingAges, bool crown, unsigned condition,
    bool branchingTimes, const Parameters& p, const Settings& settings)
{
    validate(p, settings, startAge);
    if (p.psi != 0 || p.rho != 1 || condition > 2 || branchingAges.empty())
        throw std::invalid_argument("Extant DD likelihood requires >=2 tips, psi=0, rho=1, condition 0..2");
    double last = startAge;
    for (double age : branchingAges) {
        if (!std::isfinite(age) || age <= 0 || age > last)
            throw std::invalid_argument("Invalid extant branching ages");
        last = age;
    }
    if (crown && startAge != branchingAges.front())
        throw std::invalid_argument("Crown age must equal oldest branching age");
    const std::size_t initial = crown ? 2 : 1;
    // s=0 throughout: extant identities become protected only at the present.
    std::vector<double> v(settings.maxHidden+1,0);
    v[0] = 1;
    std::size_t k = initial;
    double age = startAge, result = 0;
    for (std::size_t i = crown ? 1 : 0; i < branchingAges.size(); ++i) {
        result += interval(v,age-branchingAges[i],k,0,p,settings);
        if (p.birthModel == BirthModel::Exponential) {
            if (p.lambda0 == 0) return negativeInfinity;
            result += std::log(p.lambda0)-p.alpha*(k-1);
            for (std::size_t h=0; h<v.size(); ++h) v[h] *= std::exp(-p.alpha*h);
        } else {
            const double base = p.birth(k);
            if (base == 0) return negativeInfinity;
            result += std::log(base);
            for (std::size_t h=0; h<v.size(); ++h) v[h] *= p.birth(k+h)/base;
        }
        result += normalize(v);
        ++k;
        age = branchingAges[i];
    }
    result += interval(v,age,k,0,p,settings);
    if (v[0] <= 0) throw std::runtime_error("Extant DD endpoint underflow");
    result += std::log(v[0]);
    if (branchingTimes) result += std::lgamma(static_cast<double>(k));
    if (condition) {
        // Marked-lineage embedding: divide out the number of embeddings in a
        // surviving history. The stem weight is h+1.
        // For two crown sides the reduction gives (h+2)(h+3)/6; see derivation.
        std::fill(v.begin(),v.end(),0);
        v[0] = 1;
        const double scale = interval(v,startAge,initial,0,p,settings);
        double mass = 0;
        for (std::size_t h=0; h<v.size(); ++h) {
            if (condition == 2 && h+initial != k) continue;
            const double weight = crown ? (h+2.0)*(h+3.0)/6.0 : h+1.0;
            mass += v[h]/weight;
        }
        if (mass <= 0) throw std::runtime_error("Extant DD conditioning underflow or insufficient cap");
        result -= scale+std::log(mass);
    }
    return result;
}

// BDSTP-compatible specimen density before its labelled non-oriented shape factor.
double RevBayesCore::DDFBD::specimenLogLikelihood(double origin,
    const std::vector<SpecimenEvent>& events, std::size_t extantSamples,
    const Parameters& p, const Settings& settings)
{
    validate(p,settings,origin);
    if (p.birthModel != BirthModel::Exponential)
        throw std::invalid_argument("Specimen FBD supports exponential decline only");
    std::vector<double> v(settings.maxHidden+1,0);
    v[0]=1;
    std::size_t k=1;
    double age=origin, result=0;
    for (const auto& event : events) {
        if (!std::isfinite(event.age) || event.age<=0 || event.age>age || !k)
            throw std::invalid_argument("Invalid specimen event schedule");
        result+=interval(v,age-event.age,k,0,p,settings);
        if (event.type==SpecimenEventType::Birth) {
            if (p.lambda0==0) return negativeInfinity;
            result+=std::log(p.lambda0)-p.alpha*(k-1);
            for (std::size_t h=0;h<v.size();++h) v[h]*=std::exp(-p.alpha*h);
            ++k;
        } else {
            if (p.psi==0) return negativeInfinity;
            result+=std::log(p.psi);
            if (event.type==SpecimenEventType::TerminalSample) {
                // Sampling does not kill the lineage. Overflow is killed,
                // never reflected or counted as biological extinction.
                for (std::size_t h=v.size()-1;h>0;--h) v[h]=v[h-1];
                v[0]=0;
                --k;
            }
        }
        const double change=normalize(v);
        if (change==negativeInfinity)
            throw std::runtime_error("Specimen FBD event underflow or insufficient cutoff");
        result+=change;
        age=event.age;
    }
    if (k!=extantSamples) throw std::invalid_argument("Specimen extant count mismatch");
    result+=interval(v,age,k,0,p,settings);
    double endpoint=0, weight=1;
    for (double x:v) { endpoint+=x*weight; weight*=1-p.rho; }
    if (p.rho==0 && extantSamples) return negativeInfinity;
    if (endpoint==0) {
        // With no deaths and rho=1, any terminal fossil is impossible.
        if (p.mu==0 && p.rho==1 && std::any_of(events.begin(),events.end(),
            [](const SpecimenEvent& e){return e.type==SpecimenEventType::TerminalSample;}))
            return negativeInfinity;
        throw std::runtime_error("Specimen FBD endpoint underflow or insufficient cutoff");
    }
    return result+std::log(endpoint)+logPower(p.rho,extantSamples);
}
