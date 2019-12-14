#include "SolverFactory.h"

DBPISolver* getSolverInstance(ExperimentParams& xParams) {
    auto typeIt = dbpi_solver_types_map.find(xParams.solverID);
    if (typeIt == dbpi_solver_types_map.end()) {
        fprintf(stderr, "Invalid method type: %s\n", typeIt->second->getSolverName().c_str());
        exit(EXIT_FAILURE);
    }

    return typeIt->second->getSolverInstance(xParams);
}


template<bool naive, bool binaryAlphabet = true>
class BruteSolverFactory: public SolverFactory {
public:
    BruteSolverFactory():SolverFactory(string(naive?"naive ":"short-circuit ") + string("brute-force")) {};

    DBPISolver* getSolverInstance(ExperimentParams& xParams) {
        return new Brute_DBPI_Solver<naive, binaryAlphabet>(xParams);
    }
};

map<string, SolverFactory*> dbpi_solver_types_map =
        {{"nbf", new BruteSolverFactory<true>()},
         {"sbf", new BruteSolverFactory<false>()}
//         {"dt", "dbpit-transpone"}
         };

