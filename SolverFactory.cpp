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
class BruteSolverFactory: public SolverFactory {
public:
    BruteSolverFactory():SolverFactory(string(naive?"naive ":"short-circuit ") + string("brute-force")) {};

    DBPISolver* getSolverInstance(ExperimentParams& xParams) {
        if(xParams.isInBinaryMode())
            return new Brute_DBPI_Solver<naive, true>(xParams);
        else
            return new Brute_DBPI_Solver<naive, false>(xParams);
    }
};

map<string, SolverFactory*> dbpi_solver_types_map =
        {{NAIVE_BRUTE_FORCE_ID, new BruteSolverFactory<true>()},
         {SHORT_CIRCUIT_BRUTE_FORCE_ID, new BruteSolverFactory<false>()}
//         {"dt", "dbpit-transpone"}
         };

