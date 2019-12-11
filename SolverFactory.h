#ifndef DBPIT_SOLVERFACTORY_H
#define DBPIT_SOLVERFACTORY_H

#include "dbpi.h"
#include "xp-params.h"
#include <map>

class SolverFactory {
private:
    string solverName;

public:
    SolverFactory(string solverName): solverName(solverName) { }
    string getSolverName() { return solverName; }
    virtual DBPISolver* getSolverInstance(ExperimentParams& xParams) = 0;
};

extern map<string, SolverFactory*> dbpi_solver_types_map;

DBPISolver* getSolverInstance(ExperimentParams& xParams);

#endif //DBPIT_SOLVERFACTORY_H
