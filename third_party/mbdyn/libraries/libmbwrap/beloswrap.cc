/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
 * Paolo Mantegazza     <paolo.mantegazza@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2022(-2023) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1

  Tpetra port: AztecOO  → Belos
               Amesos   → Amesos2
               Ifpack   → Ifpack2
  This file replaces aztecoowrap.cc.
*/

#include "mbconfig.h"


#ifdef USE_TRILINOS

// include manually umfpack.h otherwisr Amesos2 is going to
// include it within extern "C"
#ifdef HAVE_UMFPACK_H
#ifdef HAVE_SUITESPARSE_EXTERN_C
#include <umfpack.h>
#else // HAVE_SUITESPARSE_EXTERN_C
extern "C" {
#include <umfpack.h>
}
#endif // HAVE_SUITESPARSE_EXTERN_C
#endif // HAVE_UMFPACK_H

#include "ls.h"
#include "linsol.h"

#ifdef USE_MPI
#include "mbcomm.h"
#endif // USE_MPI

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcpp"

#define HAVE_BLAS_SAVE HAVE_BLAS
#define HAVE_BOOL_SAVE HAVE_BOOL
#undef HAVE_BLAS
#undef HAVE_BOOL

#include "beloswrap.h"
#include "tpetravh.h"
#include "tpetraspmh.h"

/* ---- Belos (iterative solver) ----------------------------------------- */
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

/* ---- Ifpack2 (preconditioner) ----------------------------------------- */
#include <Ifpack2_Factory.hpp>

/* ---- Amesos2 (direct solver) ------------------------------------------ */
#include <Amesos2.hpp>
#include <Amesos2_Factory.hpp>

/* ---- Tpetra lifecycle ------------------------------------------------- */
#include <Tpetra_Core.hpp>
#include <Kokkos_Core.hpp>

#undef HAVE_BLAS
#undef HAVE_BOOL
#define HAVE_BLAS HAVE_BLAS_SAVE
#define HAVE_BOOL HAVE_BOOL_SAVE
#undef HAVE_BLAS_SAVE
#undef HAVE_BOOL_SAVE

#pragma GCC diagnostic pop

/* =========================================================================
 * Type shorthands
 * ========================================================================= */
using MV   = TpetraMV;
using OP   = TpetraOp;

using BelosProblem = Belos::LinearProblem<TpetraSC, MV, OP>;
using BelosSolver  = Belos::SolverManager<TpetraSC, MV, OP>;
using Ifpack2Prec  = Ifpack2::Preconditioner<TpetraSC, TpetraLO, TpetraGO, TpetraNode>;
using Amesos2Solver= Amesos2::Solver<TpetraCrs, MV>;

/* =========================================================================
 * Helper: Amesos2 solver name from flags
 * ========================================================================= */
static const char* GetAmesos2SolverName(unsigned uFlags)
{
     switch (uFlags) {
     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_UMFPACK:   return "umfpack";
     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_KLU:        return "klu2";
     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_LAPACK:     return "lapack";
     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_SUPERLU:    return "superlu";
     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_MUMPS:      return "mumps";
     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_PARDISO:    return "pardiso_mkl";
     default:
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
}

/* =========================================================================
 * Amesos2Wrapper – wraps an Amesos2::Solver with the same lazy-factorisation
 * logic the original AmesosSolver used.
 * ========================================================================= */
class Amesos2Wrapper {
public:
     Amesos2Wrapper(TpetraSparseMatrixHandler& A_,
                    TpetraVectorHandler& X_,
                    TpetraVectorHandler& B_,
                    unsigned uFlags_)
          :rA(A_), rX(X_), rB(B_),
           uFlags(uFlags_),
           bRebuildSymbolic(true),
           bRebuildNumeric(true)
     {
     }

     int Solve()
     {
          DEBUGCERR("Amesos2Wrapper::Solve()\n");

          if (pSolver.is_null()) {
               // Deferred construction: the CrsMatrix is only available
               // after the first PacMat()/EnsureFilled() cycle.
               pSolver = Amesos2::create<TpetraCrs, TpetraMV>(
                    GetAmesos2SolverName(uFlags),
                    rA.pGetTpetraCrsMatrix(),
                    rX.pGetTpetraVector(),
                    rB.pGetTpetraVector());

               if (pSolver.is_null()) {
                    silent_cerr("Trilinos/Amesos2: solver "
                                << GetAmesos2SolverName(uFlags)
                                << " is not available in this build.\n");
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
               }
          }

          do {
               if (bRebuildSymbolic) {
                    DEBUGCERR("Amesos2Wrapper::symbolicFactorization()\n");
                    pSolver->symbolicFactorization();
                    bRebuildSymbolic = false;
               }

               if (bRebuildNumeric) {
                    DEBUGCERR("Amesos2Wrapper::numericFactorization()\n");
                    try {
                         pSolver->numericFactorization();
                         bRebuildNumeric = false;
                    } catch (...) {
                         if (!bRebuildSymbolic) {
                              bRebuildSymbolic = true;
                              continue;
                         }
                         return -1;
                    }
               }
               break;
          } while (true);

          pSolver->solve();

          DEBUGCERR("Amesos2Wrapper::Solve() done\n");
          return 0;
     }

     void MatrReset()   { bRebuildNumeric  = true; }
     void MatrInitialize() { bRebuildSymbolic = bRebuildNumeric = true; pSolver.reset(); }

     TpetraSparseMatrixHandler& rA;
     TpetraVectorHandler& rX;
     TpetraVectorHandler& rB;
     unsigned uFlags;
     Teuchos::RCP<Amesos2Solver> pSolver;
     bool bRebuildSymbolic;
     bool bRebuildNumeric;
};

/* =========================================================================
 * Amesos2PrecOp — wraps an Amesos2 direct solver as a Tpetra::Operator
 * so it can serve as a right preconditioner for Belos.
 *
 * Replaces the old AztecOO + KLU/UMFPACK direct-solve-as-preconditioner
 * path: Belos sees  y = A^{-1} * x  as an operator application.
 * ========================================================================= */
class Amesos2PrecOp : public TpetraOp {
public:
     Amesos2PrecOp(const Teuchos::RCP<const TpetraMap>& pMap_,
                   const Teuchos::RCP<TpetraCrs>&       pMat_,
                   const char* solverName)
          : pMap(pMap_), pA(pMat_),
            bRebuildSymbolic(true), bRebuildNumeric(true)
     {
          pRhsMV = Teuchos::rcp(new TpetraMV(pMap, 1));
          pLhsMV = Teuchos::rcp(new TpetraMV(pMap, 1));
          pSolver = Amesos2::create<TpetraCrs, TpetraMV>(
               solverName, pA, pLhsMV, pRhsMV);
     }

     void apply(const TpetraMV& X, TpetraMV& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS,
                TpetraSC alpha = Teuchos::ScalarTraits<TpetraSC>::one(),
                TpetraSC beta  = Teuchos::ScalarTraits<TpetraSC>::zero()
                ) const override
     {
          (void)mode; // always NOTRANS for preconditioning

          // Copy input into internal RHS vector
          Tpetra::deep_copy(*pRhsMV, X);

          if (bRebuildSymbolic) {
               pSolver->symbolicFactorization();
               bRebuildSymbolic = false;
               bRebuildNumeric = true;
          }
          if (bRebuildNumeric) {
               pSolver->numericFactorization();
               bRebuildNumeric = false;
          }

          pSolver->solve();

          // Y = alpha * lhs + beta * Y
          Y.update(alpha, *pLhsMV, beta);
     }

     Teuchos::RCP<const TpetraMap> getDomainMap() const override { return pMap; }
     Teuchos::RCP<const TpetraMap> getRangeMap()  const override { return pMap; }

     void invalidateNumeric()  { bRebuildNumeric  = true; }
     void invalidateAll()      { bRebuildSymbolic = bRebuildNumeric = true; }

private:
     Teuchos::RCP<const TpetraMap> pMap;
     Teuchos::RCP<TpetraCrs>       pA;
     Teuchos::RCP<TpetraMV>        pRhsMV, pLhsMV;
     Teuchos::RCP<Amesos2Solver>   pSolver;
     mutable bool bRebuildSymbolic;
     mutable bool bRebuildNumeric;
};

/* =========================================================================
 * TpetraLinearSystem – base for both solution managers.
 *
 * Owns the communicator, the matrix handler, and both vector handlers.
 * ========================================================================= */
class TpetraLinearSystem: public SolutionManager {
public:
     explicit TpetraLinearSystem(
#ifdef USE_MPI
          MPI::Intracomm& oMpiComm,
#endif
          integer Dim);

     virtual ~TpetraLinearSystem();

#ifdef DEBUG
     virtual void IsValid() const override;
#endif

     virtual MatrixHandler* pMatHdl() const override;
     virtual VectorHandler* pResHdl() const override;
     virtual VectorHandler* pSolHdl() const override;
     virtual bool bGetConditionNumber(doublereal& dCond) const override;

protected:
     Teuchos::RCP<const TpetraComm> pComm;
     mutable TpetraVectorHandler x, b;
     mutable TpetraSparseMatrixHandler A;
};

/* -------------------------------------------------------------------------
 * TpetraLinearSystem implementation
 * ------------------------------------------------------------------------- */
TpetraLinearSystem::TpetraLinearSystem(
#ifdef USE_MPI
     MPI::Intracomm& oMpiComm,
#endif
     integer Dim)
     :
#ifdef USE_MPI
      pComm(tpetraMpiComm(oMpiComm)),
#else
      pComm(tpetraSerialComm()),
#endif
      x(Dim, pComm),
      b(Dim, pComm),
      A(Dim, Dim, 1, pComm)
{
}

TpetraLinearSystem::~TpetraLinearSystem()
{
}

#ifdef DEBUG
void TpetraLinearSystem::IsValid() const
{
     A.IsValid();
     x.IsValid();
     b.IsValid();
}
#endif

MatrixHandler* TpetraLinearSystem::pMatHdl() const { return &A; }
VectorHandler* TpetraLinearSystem::pResHdl() const { return &b; }
VectorHandler* TpetraLinearSystem::pSolHdl() const { return &x; }

bool TpetraLinearSystem::bGetConditionNumber(doublereal&) const
{
     return false;
}

/* =========================================================================
 * Amesos2SolutionManager – direct solver
 * ========================================================================= */
class Amesos2SolutionManager: public TpetraLinearSystem {
public:
     Amesos2SolutionManager(
#ifdef USE_MPI
          MPI::Intracomm& oComm,
#endif
          integer Dim,
          integer iVerbose,
          unsigned uFlags);

     virtual ~Amesos2SolutionManager();
     virtual void MatrReset()      override;
     virtual void MatrInitialize() override;
     virtual void Solve()          override;

private:
     mutable Amesos2Wrapper oSolver;
};

Amesos2SolutionManager::Amesos2SolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer Dim,
     integer iVerbose,
     unsigned uFlags)
     :TpetraLinearSystem(
#ifdef USE_MPI
          oComm,
#endif
          Dim),
      oSolver(A, x, b, uFlags)
{
}

Amesos2SolutionManager::~Amesos2SolutionManager()
{
}

void Amesos2SolutionManager::Solve()
{
     DEBUGCERR("Amesos2SolutionManager::Solve()\n");

     int ierr = oSolver.Solve();

     if (ierr != 0) {
          silent_cerr("Amesos2 error: solution failed with status " << ierr << "\n");
          throw LinearSolver::ErrFactor(-1, MBDYN_EXCEPT_ARGS);
     }
}

void Amesos2SolutionManager::MatrReset()      { oSolver.MatrReset();      }
void Amesos2SolutionManager::MatrInitialize() { oSolver.MatrInitialize(); }

/* =========================================================================
 * BelosSolutionManager – iterative solver with optional Ifpack2 preconditioner
 * ========================================================================= */
class BelosSolutionManager: public TpetraLinearSystem {
public:
     BelosSolutionManager(
#ifdef USE_MPI
          MPI::Intracomm& oComm,
#endif
          integer Dim,
          integer iMaxIter,
          doublereal dTol,
          integer iVerbose,
          unsigned uPrecondFlag);

     virtual ~BelosSolutionManager();
     virtual void Solve()     override;
     virtual void MatrReset() override;

protected:
     const integer iMaxIter;
     const doublereal dTol;
     const unsigned uPrecondFlag;

     /* Belos objects – rebuilt lazily when the matrix changes. */
     Teuchos::RCP<BelosProblem> pProblem;
     Teuchos::RCP<BelosSolver>  pSolver;
     Teuchos::RCP<Amesos2PrecOp> pAmesos2Prec;

     void BuildSolver();
     bool bSolverBuilt;
};

BelosSolutionManager::BelosSolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer Dim,
     integer iMaxIter_a,
     doublereal dTol_a,
     integer iVerbose,
     unsigned uPrecondFlag_a)
     :TpetraLinearSystem(
#ifdef USE_MPI
          oComm,
#endif
          Dim),
      iMaxIter(iMaxIter_a > 0 ? iMaxIter_a : 100),
      dTol(dTol_a > 0. ? dTol_a : 1e-10),
      uPrecondFlag(uPrecondFlag_a),
      bSolverBuilt(false)
{
}

BelosSolutionManager::~BelosSolutionManager()
{
}

void BelosSolutionManager::BuildSolver()
{
     /* Belos linear problem */
     pProblem = Teuchos::rcp(new BelosProblem(
          A.pGetTpetraCrsMatrix(),
          x.pGetTpetraVector(),
          b.pGetTpetraVector()));

     /* Ifpack2 / Amesos2 preconditioner when requested */
     if (uPrecondFlag != LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_ILUT) {
          bool bDirectSolverPrec = false;
          switch (uPrecondFlag) {
          case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_UMFPACK:
          case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_KLU:
          case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_LAPACK:
          case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_SUPERLU:
          case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_MUMPS:
          case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_PARDISO:
               bDirectSolverPrec = true;
               break;
          default:
               break;
          }

          if (bDirectSolverPrec) {
               // Use Amesos2 direct solver as preconditioner — this
               // replicates the old AztecOO + KLU/UMFPACK path where the
               // direct solver served as a preconditioner.
               pAmesos2Prec = Teuchos::rcp(new Amesos2PrecOp(
                    A.pGetTpetraCrsMatrix()->getRowMap(),
                    A.pGetTpetraCrsMatrix(),
                    GetAmesos2SolverName(uPrecondFlag)));
               pProblem->setRightPrec(pAmesos2Prec);
          } else {
               Ifpack2::Factory precFactory;
               Teuchos::RCP<Ifpack2Prec> pPrec =
                    precFactory.create("RILUK", A.pGetTpetraCrsMatrixConst());

               Teuchos::RCP<Teuchos::ParameterList> pPrecParams =
                    Teuchos::rcp(new Teuchos::ParameterList());
               pPrecParams->set("fact: iluk level-of-fill", 1);
               pPrec->setParameters(*pPrecParams);
               pPrec->initialize();
               pPrec->compute();

               pProblem->setRightPrec(pPrec);
          }
     }

     pProblem->setProblem();

     /* Belos solver parameters */
     Teuchos::RCP<Teuchos::ParameterList> pSolverParams =
          Teuchos::rcp(new Teuchos::ParameterList());
     pSolverParams->set("Maximum Iterations",    iMaxIter);
     pSolverParams->set("Convergence Tolerance", dTol);
     pSolverParams->set("Num Blocks",
                        std::max(1, std::min(iMaxIter, 300))); /* Krylov restart depth */

     Belos::SolverFactory<TpetraSC, MV, OP> factory;
     pSolver = factory.create("GMRES", pSolverParams);
     pSolver->setProblem(pProblem);

     bSolverBuilt = true;
}

void BelosSolutionManager::Solve()
{
     if (!bSolverBuilt) {
          BuildSolver();
     }

     Belos::ReturnType ret = pSolver->solve();

#ifdef DEBUG
     {
          MyVectorHandler Ax(b.iGetSize());
          A.MatVecMul(Ax, x);
          doublereal Axmb = 0., Axpb = 0.;
          for (integer i = 1; i <= b.iGetSize(); ++i) {
               Axmb += std::pow(Ax(i) - b(i), 2);
               Axpb += std::pow(Ax(i) + b(i), 2);
          }
          DEBUGCERR("||A*x - b|| / ||A*x + b|| = "
                    << std::sqrt(Axmb / Axpb) << "\n");
     }
#endif

     if (ret != Belos::Converged) {
          silent_cerr("Belos error: iterative solution did not converge\n");
          throw LinearSolver::ErrFactor(-1, MBDYN_EXCEPT_ARGS);
     }
}

void BelosSolutionManager::MatrReset()
{
     /* The matrix structure is unchanged; just invalidate the solver
      * so that the preconditioner is recomputed on the next Solve(). */
     bSolverBuilt = false;
}

/* =========================================================================
 * Public factory functions
 * ========================================================================= */
SolutionManager*
pAllocateBelosSolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer iNLD,
     integer iMaxIter,
     doublereal dTolRes,
     integer iVerbose,
     unsigned uSolverFlags)
{
     SolutionManager* pCurrSM = nullptr;
     const unsigned uPrecondFlag = uSolverFlags & LinSol::SOLVER_FLAGS_PRECOND_MASK;

#ifdef USE_MPI
     SAFENEWWITHCONSTRUCTOR(pCurrSM,
                            BelosSolutionManager,
                            BelosSolutionManager(oComm,
                                                 iNLD, iMaxIter,
                                                 dTolRes, iVerbose,
                                                 uPrecondFlag));
#else
     SAFENEWWITHCONSTRUCTOR(pCurrSM,
                            BelosSolutionManager,
                            BelosSolutionManager(iNLD, iMaxIter,
                                                 dTolRes, iVerbose,
                                                 uPrecondFlag));
#endif
     return pCurrSM;
}

SolutionManager*
pAllocateAmesos2SolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer iNLD,
     integer iVerbose,
     unsigned uSolverFlags)
{
     SolutionManager* pCurrSM = nullptr;
     const unsigned uSolverFlag = uSolverFlags & LinSol::SOLVER_FLAGS_PRECOND_MASK;

#ifdef USE_MPI
     SAFENEWWITHCONSTRUCTOR(pCurrSM,
                            Amesos2SolutionManager,
                            Amesos2SolutionManager(oComm,
                                                   iNLD, iVerbose,
                                                   uSolverFlag));
#else
     SAFENEWWITHCONSTRUCTOR(pCurrSM,
                            Amesos2SolutionManager,
                            Amesos2SolutionManager(iNLD, iVerbose,
                                                   uSolverFlag));
#endif
     return pCurrSM;
}

void
mbdyn_trilinos_finalize()
{
     // Explicitly finalize Tpetra (and Kokkos) so that Kokkos finalize
     // hooks — in particular Belos' static MultiVecPool cleanup —
     // run while all static objects are still alive.  Without this,
     // the atexit-registered Kokkos::finalize() may run after some
     // Kokkos allocation records have already been freed during C++
     // static destruction, causing SharedAllocationRecord failures.
     if (Tpetra::isInitialized()) {
          Tpetra::finalize();
     } else if (Kokkos::is_initialized()) {
          Kokkos::finalize();
     }
}

#endif  /* USE_TRILINOS */
