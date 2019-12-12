#include "SolverFactory.h"

DBPISolver* getSolverInstance(ExperimentParams& xParams) {
    auto typeIt = dbpi_solver_types_map.find(xParams.solverID);
    if (typeIt == dbpi_solver_types_map.end()) {
        fprintf(stderr, "Invalid method type: %s\n", typeIt->second->getSolverName().c_str());
        exit(EXIT_FAILURE);
    }

    return typeIt->second->getSolverInstance(xParams);
}


template<bool naive>
class NaiveBruteSolverFactory: public SolverFactory {
public:
    NaiveBruteSolverFactory():SolverFactory(string(naive?"naive ":"") + string("brute-force")) {};

    DBPISolver* getSolverInstance(ExperimentParams& xParams) {
        return new Brute_DBPI_Solver<naive>(xParams);
    }
};

map<string, SolverFactory*> dbpi_solver_types_map =
        {{"nbf", new NaiveBruteSolverFactory<true>()},
         {"sbf", new NaiveBruteSolverFactory<false>()}
//         {"dt", "dbpit-transpone"}
         };

