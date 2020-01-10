#ifndef DBPIT_SOLVERFACTORY_H
#define DBPIT_SOLVERFACTORY_H

#include "QuantizationSolver.h"
#include "xp-params.h"
#include <map>

extern map<string, SolverFactory*> dbpi_solver_types_map;

DBPISolver* getSolverInstance(ExperimentParams& xParams);

#endif //DBPIT_SOLVERFACTORY_H
