#ifndef MULTIDOMAIN_START_HH
#define MULTIDOMAIN_START_HH

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/multidomain/matrixconverter.hh>
#include <opm/multidomain/multidomainmodel.hh>
#include <opm/multidomain/multidomainproperties.hh>
#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <dune/istl/io.hh>
#include <dune/istl/solvers.hh>

namespace Opm{

template <class ProblemTypeTag>
void multidomainStart(int argc, char** argv)
{
    try {
        // Maybe initialize MPI
        Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    } catch (Dune::Exception& e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception thrown!" << std::endl;
    }

    using Simulator = GetPropType<ProblemTypeTag, Properties::Simulator>;
    using SubTypes = GetPropType<ProblemTypeTag, Properties::SubTypeTag>;
    using CouplerTypes = GetPropType<ProblemTypeTag, Properties::CouplerTypeTag>;
    using Model = GetPropType<ProblemTypeTag, Properties::Model>;
    using JacobianMatrix = typename GetPropType<ProblemTypeTag, Properties::SubTypeTag>::JacobianMatrix;
    using GlobalEqVector = typename GetPropType<ProblemTypeTag, Properties::SubTypeTag>::GlobalEqVector;
    using BlockSolutionVector = typename GetPropType<ProblemTypeTag, Properties::SubTypeTag>::GlobalSolutionVector;

    using MatrixBlock = typename GetPropType<typename SubTypes::template TypeTag<0>, Properties::SparseMatrixAdapter>::MatrixBlock::BaseType; //typename Dune::FieldMatrix<double, 2, 2>;
    using BCRSMatrix = Dune::BCRSMatrix<MatrixBlock>;
    using VectorBlock = Dune::FieldVector<double, MatrixBlock::rows>;
    using SolutionVector = Dune::BlockVector<VectorBlock>;

    // Register all parameters
    Simulator::registerParameters();

    Dune::index_constant<SubTypes::numSubDomains> numDomains;
    using namespace Dune::Hybrid;
    forEach(integralRange(numDomains), [&](const auto typeI) {
        EWOMS_END_PARAM_REGISTRATION(typename SubTypes::template TypeTag<typeI>);
    });
    Dune::index_constant<CouplerTypes::numSubCouplers> numCouplers;
    using namespace Dune::Hybrid;
    forEach(integralRange(numCouplers), [&](const auto typeI) {
        EWOMS_END_PARAM_REGISTRATION(typename CouplerTypes::template TypeTag<typeI>);
    });
    EWOMS_END_PARAM_REGISTRATION(ProblemTypeTag);

    // Initiate model
    Simulator simulator;
    auto& model = simulator.model();
    // Start time loop
    double time = 0;
    auto endTime = getPropValue<ProblemTypeTag, Properties::EndTime>();
    auto dt = getPropValue<ProblemTypeTag, Properties::InitialTimeStepSize>();
    auto maxDt = getPropValue<ProblemTypeTag, Properties::MaxTimeStepSize>();
    unsigned timeIdx = -1;
    simulator.setTimeStepSize(dt);
    simulator.setTime(-dt, timeIdx);
    model.applyInitialSolution();
    model.writeOutput();
    simulator.setTime(0.0, ++timeIdx);

    double simulationTime = 0;
    double solverTimeTot = 0;
    double linearizationTimeTot = 0;
    int maxNewtonIt = 25;
    while (time <= endTime - dt + 1e-8) {
        std::cout << "Starting time step\n"
                  << "t = " << time << ", dt = " << dt << '\n';
        auto start_time = std::chrono::high_resolution_clock::now();
        simulator.setTime(time, timeIdx);
        BlockSolutionVector& nextSolution = model.solution(/*historyIdx=*/0);
        BlockSolutionVector currentSolution(nextSolution);
        BlockSolutionVector oldSolution(currentSolution);

        double error = 1e10;

        double solverTime = 0;
        double linearizationTime = 0;
        int newtonIterationIdx = 0;

        forEach(integralRange(numDomains), [&](const auto domainI) {
            model.template model<domainI>().newtonMethod().setIterationIndex(newtonIterationIdx);
        });

        while (error > 1e-8 && newtonIterationIdx < maxNewtonIt) {
            // BlockSolutionVector& oldSolution = model.solution(/*historyIdx=*/1);
            // auto error2 = oldSolution.two_norm();
            // std::cout << "solution norm: " << error2<<std::endl;

            // make the current solution to the old one
            auto start_lin = std::chrono::high_resolution_clock::now();
            currentSolution = nextSolution;
            try {
                model.linearizer().linearize();
            } catch (Opm::NumericalIssue) {
                newtonIterationIdx = maxNewtonIt;
                break;
            }
            auto end_lin = std::chrono::high_resolution_clock::now();
            linearizationTime += std::chrono::duration_cast<std::chrono::milliseconds>(end_lin - start_lin).count();
            auto& jac = model.linearizer().jacobian();
            auto& res = model.linearizer().residual();
            GlobalEqVector blockUpdate(res);
            // Set up linear solver. First convert jacobian and residual to simple BCRS matrix
            // and block vector.
            const auto& bcrs = Opm::MatrixConverter<JacobianMatrix, MatrixBlock>::multiTypeToBCRSMatrix(jac);
            const auto& blockRes = Opm::VectorConverter<GlobalEqVector, VectorBlock>::multiTypeToBlockVector(res);

            double norm = 0.0;
            Dune::index_constant<0> _0;
            for (const auto& row : jac[_0][_0])
                for (const auto& entry : row)
                    norm += entry.infinity_norm();

            error = blockRes.infinity_norm();
            if (error < 1e-9 && newtonIterationIdx > 0) {
                std::cout << "Newton iteration finished. Error: " << error << std::endl;
                break;
            }

            auto start_solve = std::chrono::high_resolution_clock::now();

            SolutionVector solutionUpdate(blockRes.size());
            SolutionVector residual(blockRes);

            using LinearOperator = Dune::MatrixAdapter<BCRSMatrix, SolutionVector, SolutionVector>;
            using CoarsenCriterion = typename Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<BCRSMatrix, Dune::Amg::FrobeniusNorm>>;
            ///typedef typename Dune::Amg::SymmetricCriterion<BCRSMatrix, Dune::Amg::FrobeniusNorm> CCC;
            using ILU0 = Dune::SeqILU0<BCRSMatrix, SolutionVector, SolutionVector>;
            using SOR = Dune::SeqSOR<BCRSMatrix, SolutionVector, SolutionVector>;
            using AMG = typename Dune::Amg::AMG<LinearOperator, SolutionVector, SOR>;
            using SmootherArgs = typename Dune::Amg::SmootherTraits<SOR>::Arguments;

            Dune::InverseOperatorResult statistics;
            CoarsenCriterion coarsenCriterion(/*maxLevel=*/15, 5000);
            coarsenCriterion.setMinCoarsenRate(1.05);
            coarsenCriterion.setDebugLevel(0); // make the AMG shut up
            // coarsenCriterion.setAccumulate(Dune::Amg::noAccu);
            coarsenCriterion.setAccumulate(Dune::Amg::atOnceAccu);
            coarsenCriterion.setSkipIsolated(false);
            SmootherArgs smootherArgs;
            smootherArgs.iterations = 1;
            smootherArgs.relaxationFactor = 1.0;

            LinearOperator linearOperator(bcrs);
            ILU0 preconditioner(bcrs, 1.0);
            // AMG preconditioner(linearOperator, coarsenCriterion, smootherArgs);
            // // Define solver and solve
            //Dune::UMFPack<BCRSMatrix> solver(bcrs);
            Dune::BiCGSTABSolver<SolutionVector> solver(linearOperator, preconditioner, 1e-4, 400, 0); //-25,

            solutionUpdate = 0.0;
            try {
                solver.apply(solutionUpdate, residual, statistics);
            } catch (...) {
                newtonIterationIdx = maxNewtonIt;
                break;
            }
            auto end_solve = std::chrono::high_resolution_clock::now();
            solverTime += std::chrono::duration_cast<std::chrono::milliseconds>(end_solve - start_solve).count();
            // Format solution back to a multitype block vector
            blockUpdate = 0.0;
            Opm::VectorConverter<GlobalEqVector, VectorBlock>::retrieveValues(blockUpdate, solutionUpdate);

            // Update newton iteration index. Must be done before we possible do the line search.
            newtonIterationIdx++;
            forEach(integralRange(numDomains), [&](const auto typeI) {
                model.template model<typeI>().newtonMethod().setIterationIndex(newtonIterationIdx);
            });
            double beta = 1.0;
            double error_search;
            for (int i = 0; i < 3; i++) { // line search

                // Update solution
                using namespace Dune::Hybrid;
                forEach(integralRange(Dune::Hybrid::size(nextSolution)), [&](const auto domainI) {
                    // update the DOFs of the auxiliary equations
                    size_t numDof = model.template model<domainI>().numTotalDof();
                    for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
                        for (size_t phaseIdx = 0; phaseIdx < 2; ++phaseIdx) {
                            nextSolution[domainI][dofIdx][phaseIdx] = currentSolution[domainI][dofIdx][phaseIdx];
                            nextSolution[domainI][dofIdx][phaseIdx] -= beta * blockUpdate[domainI][dofIdx][phaseIdx];
                        }
                    }
                });

                model.setSolution(nextSolution);

                // make sure that the intensive quantities get recalculated at the next
                // linearization
                model.invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);

                try {
                    model.linearizer().linearize();
                } catch (Opm::NumericalIssue) {
                    newtonIterationIdx = maxNewtonIt;
                    break;
                }

                res = model.linearizer().residual();
                error_search = res.infinity_norm();
                if (i > 0)
                    std::cout << "Reducing step length: beta = " << beta << ", error = " << error_search << std::endl;
                if (error_search < error) {
                    break;
                }
                beta /= 3.0;
            }

            error = error_search;
            std::cout << "Newton iteration finished. Error: " << error << std::endl;
            if (std::isnan(error) || error > 1e8) {
                newtonIterationIdx = maxNewtonIt;
                break;
            }
        }
        if (newtonIterationIdx >= maxNewtonIt) {
            dt *= 0.5;
            simulator.setTimeStepSize(dt);
            model.setSolution(oldSolution);
            model.invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
            std::cout << "Failed Newton. reducing time step.\n";
            continue;
        }

        linearizationTimeTot += linearizationTime;
        solverTimeTot += solverTime;
        std::cout << "Time step took: " << (linearizationTime + solverTime) / 1000.0 << '\n'
                  << "   Linearization: " << linearizationTime / 1000.0 << '\n'
                  << "   Solve: " << solverTime / 1000.0 << '\n';

        std::cout << "write output" << std::endl;
        model.writeOutput();
        std::cout << "Advancing time step." << std::endl;
        model.advanceTimeLevel();
        time += dt;
        timeIdx++;
        if (dt < maxDt - 1e-5 && newtonIterationIdx < 7)
            dt *= 1.2;
        simulator.setTimeStepSize(dt);
        std::cout << "Finished time step." << std::endl;

        auto end_time = std::chrono::high_resolution_clock::now();
        simulationTime += std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    }
    std::cout << "Simulation took " << simulationTime / 1000.0 << " s\n"
              << "   Linearization: " << linearizationTimeTot / 1000.0 << " s\n"
              << "   Solve: " << solverTimeTot / 1000.0 << " s\n";
}

} // end namespace Opm

#endif // MULTIDOMAIN_START_HH