// Standalone numerical regression driver; links only the DDFBD kernel.
#include "DiversityDependentFbdLikelihood.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace RevBayesCore::DDFBD;
int main(int argc, char** argv) {
    try {
        if (argc < 10) throw std::invalid_argument("insufficient arguments");
        const std::string mode = argv[1];
        Parameters p{std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]),
                     std::stod(argv[5]), std::stod(argv[6])};
        const double origin = std::stod(argv[7]);
        Settings settings;
        settings.maxHidden = std::stoull(argv[8]);
        settings.tolerance = std::stod(argv[9]);
        double answer;
        if (mode == "none" || mode == "noextant" || mode == "logsampling" || mode == "logextant") {
            if (argc != 10) throw std::invalid_argument("extra conditioning arguments");
            answer = (mode == "logsampling" || mode == "logextant")
                ? logObservationProbability(origin, settings.maxHidden, mode == "logsampling", p, settings)
                : noObservationProbability(origin, settings.maxHidden, mode == "none", p, settings);
        } else if (mode == "ll") {
            if (argc < 14) throw std::invalid_argument("missing tree arguments");
            const auto fossils = std::stoull(argv[10]);
            const auto sampled = std::stoull(argv[11]);
            const auto unsampled = std::stoull(argv[12]);
            const auto count = std::stoull(argv[13]);
            if (static_cast<std::size_t>(argc) != 14+2*count) throw std::invalid_argument("event argument count mismatch");
            std::vector<Event> events;
            for (std::size_t i=0; i<count; ++i) {
                const std::string type = argv[15+2*i];
                EventType eventType;
                if (type == "B") eventType=EventType::Birth;
                else if (type == "P") eventType=EventType::Protect;
                else if (type == "D") eventType=EventType::Death;
                else throw std::invalid_argument("unknown event type");
                events.push_back({std::stod(argv[14+2*i]), eventType});
            }
            answer=logLikelihood(origin, events, fossils, sampled, unsampled, p, settings);
        } else throw std::invalid_argument("unknown mode");
        std::cout << std::setprecision(17) << answer << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
