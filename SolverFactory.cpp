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
        else {
            switch(xParams.bytesPerElement) {
                case 1: return new Brute_DBPI_Solver<naive, false, uint8_t>(xParams);
                case 2: return new Brute_DBPI_Solver<naive, false, uint16_t>(xParams);
                case 4: return new Brute_DBPI_Solver<naive, false, uint32_t>(xParams);
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        }
    }
};

class QuantizationBasedSolverFactory: public SolverFactory {
public:
    QuantizationBasedSolverFactory():SolverFactory(string("quantization-based")) {};

    DBPISolver* getSolverInstance(ExperimentParams& xParams) {
        if(xParams.isInBinaryMode()) {
            fprintf(stderr, "Quantization unsupported for binary dataset.\n");
            exit(EXIT_FAILURE);
        } else {
            DBPISolver* postSolver = BruteSolverFactory<false>().getSolverInstance(xParams);
            SolverFactory* binarySolverFactory = new BruteSolverFactory<false>();
            switch(xParams.bytesPerElement) {
                case 1: return new QuantizationBased_DBPI_Solver<uint8_t>(xParams, new SimpleBinaryQuantizer<uint8_t>(), binarySolverFactory, postSolver);
                case 2: return new QuantizationBased_DBPI_Solver<uint16_t>(xParams, new SimpleBinaryQuantizer<uint16_t>(), binarySolverFactory, postSolver);
                case 4: return new QuantizationBased_DBPI_Solver<uint32_t>(xParams, new SimpleBinaryQuantizer<uint32_t>(), binarySolverFactory, postSolver);
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        }
    }
};


map<string, SolverFactory*> dbpi_solver_types_map =
        {{NAIVE_BRUTE_FORCE_ID, new BruteSolverFactory<true>()},
         {SHORT_CIRCUIT_BRUTE_FORCE_ID, new BruteSolverFactory<false>()},
         {QUATIZATION_BASED_SOLVER_ID, new QuantizationBasedSolverFactory()}
         };

