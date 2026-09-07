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

  Tpetra port:
    NOX::Epetra::*              -> NOX::Thyra::*
    AztecOO                     -> Belos (via Stratimikos)
    Epetra_MpiComm/SerialComm   -> Tpetra communicator
    Epetra_Map/Operator/Vector  -> Thyra/Tpetra equivalents

  The linearisation interface (computeF / computeJacobian /
  computePreconditioner) is re-expressed through the NOX::Thyra model
  evaluator thin wrapper ModelEvaluatorWrapper defined below.
*/

#include "mbconfig.h"

#ifdef USE_TRILINOS
#include <set>
#include "solman.h"
#include "solver.h"
#include "noxsolver.h"
#include "output.h"
// Note: tpetraspmh.h is NOT included here intentionally.
// ModelEvaluatorWrapper::create_W_op() uses MBDynJacobianOp which delegates
// to the abstract MatrixHandler interface, so no TpetraSparseMatrixHandler
// is needed in the non-matrix-free path.
#ifdef USE_MPI
#include "mbcomm.h"
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcpp"

#define HAVE_BLAS_SAVE HAVE_BLAS
#define HAVE_BOOL_SAVE HAVE_BOOL
#undef HAVE_BLAS
#undef HAVE_BOOL

/* ---------- Tpetra / Thyra core ---------------------------------------- */
#include "tpetra_types.h"
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryBase.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpSourceBase.hpp>
#include <Thyra_DefaultPreconditioner.hpp>

/* ---------- Stratimikos (Belos via Thyra) ------------------------------- */
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>

/* ---------- NOX Thyra interface ----------------------------------------- */
#include <NOX.H>
#include <NOX_Thyra.H>
#include <NOX_Thyra_Group.H>
#include <NOX_Thyra_Vector.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Solver_LineSearchBased.H>

/* ---------- Thyra model evaluator --------------------------------------- */
#include <Thyra_ModelEvaluator.hpp>
#include <Thyra_StateFuncModelEvaluatorBase.hpp>

/* ---------- Teuchos ----------------------------------------------------- */
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <BelosOutputManager.hpp>
#include <BelosSolverFactory.hpp>

#undef HAVE_BLAS
#undef HAVE_BOOL
#define HAVE_BLAS HAVE_BLAS_SAVE
#define HAVE_BOOL HAVE_BOOL_SAVE
#undef HAVE_BLAS_SAVE
#undef HAVE_BOOL_SAVE

#pragma GCC diagnostic pop

#ifdef DEBUG_JACOBIAN
#include "sp_gradient_spmh.h"
#endif

/* -------------------------------------------------------------------------
 * NoxSolverParameters
 * ------------------------------------------------------------------------- */
NoxSolverParameters::NoxSolverParameters()
     :CommonNonlinearSolverParam(SOLVER_LINESEARCH_BASED   |
                                 JACOBIAN_NEWTON           |
                                 DIRECTION_NEWTON          |
                                 FORCING_TERM_CONSTANT     |
                                 LINESEARCH_BACKTRACK      |
                                 LINEAR_SOLVER_BLOCK_GMRES |
                                 RECOVERY_STEP_TYPE_CONST,
                                 0,
                                 false),
      dWrmsRelTol(0.),
      dWrmsAbsTol(0.),
      dTolLinSol(1e-10),
      dMinStep(1e-12),
      dRecoveryStep(1.),
      dForcingTermMinTol(1e-6),
      dForcingTermMaxTol(1e-2),
      dForcingTermAlpha(1.5),
      dForcingTermGamma(0.9),
      iMaxIterLinSol(1000),
      iKrylovSubSpaceSize(300),
      iMaxIterLineSearch(200),
      iInnerIterBeforeAssembly(std::numeric_limits<integer>::max())
{
}

/* -------------------------------------------------------------------------
 * eGetBelosParamType – query the type of a Belos solver parameter.
 *
 * Creates a solver instance with null parameters the first time each
 * solver-type name is requested, caches its getValidParameters() list, and
 * returns a BelosParamType enum value that the input parser can use to
 * dispatch to the correct HighParser::Get*() method.
 *
 * Thread-safety note: the cache is a function-local static protected by the
 * C++11 guarantee of initialisation on first use; subsequent accesses are
 * read-only and require no locking.
 * ------------------------------------------------------------------------- */
BelosParamType
eGetBelosParamType(const std::string& sSolverType, const std::string& sParamName)
{
     // Cache: solver-type name → its full valid-parameter list.
     // Built lazily; each entry is created exactly once.
     using ValidParamCache =
          std::map<std::string,
                   Teuchos::RCP<const Teuchos::ParameterList>>;
     static ValidParamCache oCache;

     auto it = oCache.find(sSolverType);
     if (it == oCache.end()) {
          // Create a throw-away solver just to interrogate its valid params.
          //
          // In some Trilinos builds Belos::SolverFactory<SC,MV,OP> is a
          // typedef alias for Belos::Impl::SolverFactoryParent<SC,MV,OP>,
          // whose default constructor is protected.  We work around this by
          // deriving a minimal local subclass: a derived-class constructor is
          // permitted to call the protected base-class constructor.
          struct TBelosFactory
               : public Belos::SolverFactory<TpetraSC, TpetraMV, TpetraOp> {
               TBelosFactory() = default;
          };
          TBelosFactory factory;
          Teuchos::RCP<Belos::SolverManager<TpetraSC, TpetraMV, TpetraOp>> pSolver;
          try {
               pSolver = factory.create(sSolverType,
                                        Teuchos::null /* use defaults */);
          } catch (...) {
               // Unknown solver type – cannot proceed.
               return BelosParamType::UNKNOWN;
          }
          it = oCache.emplace(sSolverType, pSolver->getValidParameters()).first;
     }

     const Teuchos::ParameterList& oValid = *it->second;

     if (!oValid.isParameter(sParamName)) {
          return BelosParamType::UNKNOWN;
     }

     const Teuchos::ParameterEntry& oEntry = oValid.getEntry(sParamName);

     // Dispatch on the four scalar types that Belos actually uses for
     // user-settable parameters.  Anything else (RCP<ostream>, Array<…>,
     // etc.) is flagged as UNSETTABLE so the parser can emit a clear error.
     if (oEntry.isType<bool>())        return BelosParamType::BOOL;
     if (oEntry.isType<int>())         return BelosParamType::INT;
     if (oEntry.isType<double>())      return BelosParamType::DOUBLE;
     if (oEntry.isType<std::string>()) return BelosParamType::STRING;

     return BelosParamType::UNSETTABLE;
}

namespace {

class NoxNonlinearSolver;

/* =========================================================================
 * Status tests
 * ========================================================================= */
class NoxStatusTest : public NOX::StatusTest::Generic {
public:
     explicit NoxStatusTest(NoxNonlinearSolver& s)
          : oNoxSolver(s), eStatus(NOX::StatusTest::Unevaluated) {}
     virtual ~NoxStatusTest() {}
     NOX::StatusTest::StatusType getStatus() const override { return eStatus; }
     void Reset() { eStatus = NOX::StatusTest::Unevaluated; }
protected:
     NoxNonlinearSolver& oNoxSolver;
     NOX::StatusTest::StatusType eStatus;
};

class NoxResidualTest : public NoxStatusTest {
public:
     explicit NoxResidualTest(NoxNonlinearSolver& s)
          : NoxStatusTest(s), dErrRes(-1.), dErrResDiff(-1.), dTolRes(-2.) {}
     ~NoxResidualTest() {}
     NOX::StatusTest::StatusType
     checkStatus(const NOX::Solver::Generic& problem,
                 NOX::StatusTest::CheckType checkType) override;
     std::ostream& print(std::ostream& stream, int indent) const override;
     void Reset() { NoxStatusTest::Reset(); dErrRes = dErrResDiff = -1.; }
     void SetTolerance(doublereal dTol) { ASSERT(dTol >= 0.); dTolRes = dTol; }
     doublereal dGetTest()     const { return dErrRes;     }
     doublereal dGetTestDiff() const { return dErrResDiff; }
private:
     doublereal dErrRes, dErrResDiff, dTolRes;
};

class NoxSolutionTest : public NoxStatusTest {
public:
     explicit NoxSolutionTest(NoxNonlinearSolver& s)
          : NoxStatusTest(s), dErrSol(-1.), dTolSol(-2.) {}
     ~NoxSolutionTest() {}
     NOX::StatusTest::StatusType
     checkStatus(const NOX::Solver::Generic& problem,
                 NOX::StatusTest::CheckType checkType) override;
     std::ostream& print(std::ostream& stream, int indent) const override;
     void Reset() { NoxStatusTest::Reset(); dErrSol = -1.; }
     void SetTolerance(doublereal dTol) { ASSERT(dTol >= 0.); dTolSol = dTol; }
     doublereal dGetTolerance() const { ASSERT(dTolSol >= 0.); return dTolSol; }
     doublereal dGetTest()      const { return dErrSol; }
private:
     doublereal dErrSol, dTolSol;
};

/* =========================================================================
 * TpetraMatFreeJacOper
 *
 * Matrix-free Jacobian operator for the JACOBIAN_NEWTON_KRYLOV path.
 * Replaces NoxMatrixFreeJacOper (Epetra_Operator) from the original.
 * Implements Thyra::LinearOpBase<SC> so it is handed directly to the
 * NOX::Thyra::Group as the Jacobian and applies J*x on the fly via
 * pNonlinearProblem->Jacobian(&Y, &X).
 * ========================================================================= */
class TpetraMatFreeJacOper : public Thyra::LinearOpBase<TpetraSC>
{
public:
     using SC = TpetraSC;

     TpetraMatFreeJacOper(
          NoxNonlinearSolver& solver_a,
          const Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>& pSpace_a)
          : oNoxSolver(solver_a), pSpace(pSpace_a)
#ifdef DEBUG_JACOBIAN
          , pA(nullptr)
#endif
     {}

     ~TpetraMatFreeJacOper()
     {
#ifdef DEBUG_JACOBIAN
          if (pA) { SAFEDELETE(pA); }
#endif
     }

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     range()  const override { return pSpace; }

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     domain() const override { return pSpace; }

     bool opSupportedImpl(Thyra::EOpTransp M_trans) const override
     {
          return (M_trans == Thyra::NOTRANS);
     }

protected:
     void applyImpl(
          const Thyra::EOpTransp                          M_trans,
          const Thyra::MultiVectorBase<SC>&               X_in,
          const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
          const SC                                        alpha,
          const SC                                        beta) const override;

private:
     NoxNonlinearSolver& oNoxSolver;
     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>> pSpace;
#ifdef DEBUG_JACOBIAN
     mutable SpGradientSparseMatrixHandler* pA;
     mutable MyVectorHandler AX;
#endif
};

/* =========================================================================
 * MBDynJacobianOp
 *
 * A Thyra::LinearOpBase<SC> that applies the Jacobian stored in MBDyn's
 * MatrixHandler (via MatVecMul / MatTVecMul).  Works with *any* matrix
 * type that MBDyn's SolutionManager may wrap — no assumption of Tpetra.
 *
 * This replaces the Tpetra-specific Thyra::TpetraLinearOp path that
 * required dynamic_cast<TpetraSparseMatrixHandler*>.
 * ========================================================================= */
class MBDynJacobianOp : public Thyra::LinearOpBase<TpetraSC>
{
public:
     using SC = TpetraSC;

     MBDynJacobianOp(NoxNonlinearSolver& solver_a,
                     const Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>& pSpace_a)
          : oNoxSolver(solver_a), pSpace(pSpace_a) {}

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     range()  const override { return pSpace; }

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     domain() const override { return pSpace; }

     bool opSupportedImpl(Thyra::EOpTransp M_trans) const override
     {
          return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS);
     }

protected:
     void applyImpl(
          const Thyra::EOpTransp                          M_trans,
          const Thyra::MultiVectorBase<SC>&               X_in,
          const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
          const SC                                        alpha,
          const SC                                        beta) const override;

private:
     NoxNonlinearSolver& oNoxSolver;
     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>> pSpace;
};

/* =========================================================================
 * MBDynLinearOpWithSolve
 *
 * A Thyra::LinearOpWithSolveBase<SC> that:
 *   - applies  J   via MBDynJacobianOp::applyImpl  (MatVecMul)
 *   - solves   J x = b via pSolutionManager->Solve()
 *
 * This is the Tpetra-port equivalent of the Epetra ApplyInverse path.
 * It works with *any* SolutionManager / linear solver that MBDyn supports.
 * ========================================================================= */
class MBDynLinearOpWithSolve : public Thyra::LinearOpWithSolveBase<TpetraSC>
{
public:
     using SC = TpetraSC;

     MBDynLinearOpWithSolve(NoxNonlinearSolver& solver_a,
                            const Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>& pSpace_a)
          : oNoxSolver(solver_a), pSpace(pSpace_a), pJacOp(Teuchos::null) {}

     /* ---- LinearOpBase interface ---------------------------------------- */
     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     range()  const override { return pSpace; }

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     domain() const override { return pSpace; }

     bool opSupportedImpl(Thyra::EOpTransp M_trans) const override
     {
          return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS);
     }

     /* ---- LinearOpWithSolveBase interface -------------------------------- */
     bool solveSupportsImpl(Thyra::EOpTransp M_trans) const override
     {
          return (M_trans == Thyra::NOTRANS);
     }

     bool solveSupportsSolveMeasureTypeImpl(
          Thyra::EOpTransp                      /* M_trans */,
          const Thyra::SolveMeasureType& /* solveMeasureType */) const override
     {
          return true;
     }

protected:
     /* ---- apply: delegate to MBDyn MatrixHandler (MatVecMul) ------------ */
     void applyImpl(
          const Thyra::EOpTransp                          M_trans,
          const Thyra::MultiVectorBase<SC>&               X_in,
          const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
          const SC                                        alpha,
          const SC                                        beta) const override;

     /* ---- solveImpl: delegate to MBDyn SolutionManager->Solve() --------- */
     Thyra::SolveStatus<SC> solveImpl(
          const Thyra::EOpTransp                          M_trans,
          const Thyra::MultiVectorBase<SC>&               B,
          const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& X,
          Teuchos::Ptr<const Thyra::SolveCriteria<SC>>    solveCriteria)
          const override;

private:
     NoxNonlinearSolver& oNoxSolver;
     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>> pSpace;
     mutable Teuchos::RCP<Thyra::LinearOpBase<SC>>  pJacOp; // unused; kept for symmetry
};

/* =========================================================================
 * MBDynLOWSFactory
 *
 * A Thyra::LinearOpWithSolveFactoryBase<SC> that creates
 * MBDynLinearOpWithSolve objects.  NOX::Thyra::Group calls
 * createOp() / initializeOp() to get the object it will use for
 * linear solves; by returning an MBDynLinearOpWithSolve we bypass
 * Stratimikos/Belos/Ifpack2 entirely for the non-matrix-free path,
 * and instead delegate directly to pSolutionManager->Solve() — exactly
 * what the Epetra ApplyInverse path did.
 * ========================================================================= */
class MBDynLOWSFactory
     : public Thyra::LinearOpWithSolveFactoryBase<TpetraSC>
{
public:
     using SC = TpetraSC;

     MBDynLOWSFactory(NoxNonlinearSolver& solver_a,
                      const Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>& pSpace_a)
          : oNoxSolver(solver_a), pSpace(pSpace_a) {}

     bool isCompatible(const Thyra::LinearOpSourceBase<SC>& /* fwdOpSrc */) const override
     { return true; }

     Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC>> createOp() const override
     {
          return Teuchos::rcp(new MBDynLinearOpWithSolve(oNoxSolver, pSpace));
     }

     void initializeOp(
          const Teuchos::RCP<const Thyra::LinearOpSourceBase<SC>>& /* fwdOpSrc */,
          Thyra::LinearOpWithSolveBase<SC>*                          Op,
          const Thyra::ESupportSolveUse                            /* supportSolveUse */)
          const override
     {
          // Nothing to do: MBDynLinearOpWithSolve always reads the current
          // pSolutionManager state at solve time.
          (void)Op;
     }

     void uninitializeOp(
          Thyra::LinearOpWithSolveBase<SC>*                          /* Op */,
          Teuchos::RCP<const Thyra::LinearOpSourceBase<SC>>*         /* fwdOpSrc */,
          Teuchos::RCP<const Thyra::PreconditionerBase<SC>>*         /* prec */,
          Teuchos::RCP<const Thyra::LinearOpSourceBase<SC>>*         /* approxFwdOpSrc */,
          Thyra::ESupportSolveUse*                                   /* supportSolveUse */)
          const override
     {}

     bool supportsPreconditionerInputType(
          const Thyra::EPreconditionerInputType /* precOpType */) const override
     { return true; }

     void initializePreconditionedOp(
          const Teuchos::RCP<const Thyra::LinearOpSourceBase<SC>>&    /* fwdOpSrc */,
          const Teuchos::RCP<const Thyra::PreconditionerBase<SC>>&    /* prec */,
          Thyra::LinearOpWithSolveBase<SC>*                            /* Op */,
          const Thyra::ESupportSolveUse                               /* supportSolveUse */)
          const override
     {
          // Nothing to do: MBDynLinearOpWithSolve always reads the current
          // pSolutionManager state at solve time, and the preconditioner
          // (MBDynPrecOp) similarly delegates to pSolutionManager->Solve().
     }

     void initializeApproxPreconditionedOp(
          const Teuchos::RCP<const Thyra::LinearOpSourceBase<SC>>&    /* fwdOpSrc */,
          const Teuchos::RCP<const Thyra::LinearOpSourceBase<SC>>&    /* approxFwdOpSrc */,
          Thyra::LinearOpWithSolveBase<SC>*                            /* Op */,
          const Thyra::ESupportSolveUse                               /* supportSolveUse */)
          const override
     {
          // Same as above.
     }

     std::string description() const override { return "MBDynLOWSFactory"; }

     /* Teuchos::ParameterListAcceptor stubs (required interface) */
     void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>&) override {}
     Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList() override
     { return Teuchos::null; }
     Teuchos::RCP<Teuchos::ParameterList>       unsetParameterList() override
     { return Teuchos::null; }
     Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const override
     { return Teuchos::null; }
     Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override
     { return Teuchos::null; }

private:
     NoxNonlinearSolver& oNoxSolver;
     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>> pSpace;
};

/* =========================================================================
 * MBDynPrecOp
 *
 * Preconditioner operator for the JFNK path.  Its apply() computes
 * y = J^{-1} * x  by delegating to pSolutionManager->Solve(), exactly as
 * NoxNonlinearSolver::ApplyInverse did in the Epetra version.
 *
 * Belos calls apply() on the preconditioner at every GMRES iteration.
 * With the assembled Jacobian as an exact (or near-exact) preconditioner,
 * GMRES converges in very few iterations.
 * ========================================================================= */
class MBDynPrecOp : public Thyra::LinearOpBase<TpetraSC>
{
public:
     using SC = TpetraSC;

     MBDynPrecOp(NoxNonlinearSolver& solver_a,
                 const Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>& pSpace_a)
          : oNoxSolver(solver_a), pSpace(pSpace_a) {}

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     range()  const override { return pSpace; }

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     domain() const override { return pSpace; }

protected:
     bool opSupportedImpl(Thyra::EOpTransp M_trans) const override
     {
          return (M_trans == Thyra::NOTRANS);
     }

     void applyImpl(
          const Thyra::EOpTransp                          M_trans,
          const Thyra::MultiVectorBase<SC>&               X_in,
          const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
          const SC                                        alpha,
          const SC                                        beta) const override;

private:
     friend class NoxNonlinearSolver;

     NoxNonlinearSolver& oNoxSolver;
     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>> pSpace;
};

/* =========================================================================
 * ModelEvaluatorWrapper
 *
 * Bridges MBDyn residual / Jacobian to NOX::Thyra via Thyra::ModelEvaluator.
 *
 * Key design points vs. the Epetra version:
 *
 *  1. create_W_op() returns a MBDynJacobianOp that delegates
 *     matrix-vector products to pSolutionManager->pMatHdl() via the
 *     abstract MatrixHandler interface — no assumption of TpetraSparseMatrix.
 *
 *  2. For JACOBIAN_NEWTON_KRYLOV, create_W_op() returns a
 *     TpetraMatFreeJacOper instead.
 *
 *  3. evalModelImpl fully mirrors computeF + computeJacobian from the
 *     original, including:
 *       - residual branch         <- computeF
 *       - W_op branch             <- computeJacobian / Jacobian()
 *       - W_prec branch           <- computePreconditioner / Jacobian()
 *       - MBDyn convention that Residual() must precede Jacobian()
 *       - line-search lambda and iteration-counter bookkeeping
 *       - NOX sign convention (residual negated)
 *
 *  4. The linear solve is performed by MBDynLinearOpWithSolve, which
 *     delegates to pSolutionManager->Solve() — exactly as ApplyInverse
 *     did in the Epetra version, and works with any SolutionManager.
 *
 *  5. For JACOBIAN_NEWTON_KRYLOV, create_W_prec() returns a
 *     Thyra::DefaultPreconditioner wrapping MBDynPrecOp, and
 *     evalModelImpl assembles the Jacobian matrix when W_prec is
 *     requested — mirroring the Epetra computePreconditioner callback.
 * ========================================================================= */
class ModelEvaluatorWrapper
     : public Thyra::StateFuncModelEvaluatorBase<TpetraSC>
{
public:
     using SC = TpetraSC;

     ModelEvaluatorWrapper(NoxNonlinearSolver& solver,
                           integer iSize,
                           const Teuchos::RCP<const TpetraComm>& pComm_a);

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     get_x_space() const override { return pSpace; }

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     get_f_space() const override { return pSpace; }

     Thyra::ModelEvaluatorBase::InArgs<SC>
     getNominalValues() const override { return oNominalValues; }

     Thyra::ModelEvaluatorBase::InArgs<SC>  createInArgs()      const override;
     Thyra::ModelEvaluatorBase::OutArgs<SC> createOutArgsImpl()  const override;

     Teuchos::RCP<Thyra::LinearOpBase<SC>> create_W_op() const override;
     Teuchos::RCP<Thyra::PreconditionerBase<SC>> create_W_prec() const override;

     void evalModelImpl(
          const Thyra::ModelEvaluatorBase::InArgs<SC>&  inArgs,
          const Thyra::ModelEvaluatorBase::OutArgs<SC>& outArgs) const override;

     void Rebuild(integer iSize);

     const Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>&
     GetSpace() const { return pSpace; }

     // Returns the LOWS factory for NOX::Thyra::Group.
     // For the non-matrix-free path this is MBDynLOWSFactory.
     // For the matrix-free path a Stratimikos factory is still used.
     Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC>>
     GetLOWSFactory() const { return pLOWSFactory; }

private:
     friend class NoxNonlinearSolver;

     NoxNonlinearSolver& oNoxSolver;
     Teuchos::RCP<const TpetraComm>                             pComm;
     Teuchos::RCP<const TpetraMap>                              pMap;
     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>             pSpace;
     Thyra::ModelEvaluatorBase::InArgs<SC>                      oNominalValues;
     mutable Teuchos::RCP<Thyra::LinearOpBase<SC>>              pJacobianOp;
     mutable Teuchos::RCP<Thyra::PreconditionerBase<SC>>        pPreconditioner;
     Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC>>      pLOWSFactory;
};

/* =========================================================================
 * MBDynThyraGroup — work around NOX::Thyra::Group shared-Jacobian
 *                   ownership bug in operator= and clone(DeepCopy)
 *
 * NOX::Thyra::Group::operator= checks  this->isJacobian()  which
 * requires shared_jacobian_->isOwner(this) — but ownership has not
 * been transferred yet at that point.  The Epetra equivalent
 * (NOX::Epetra::Group) correctly checks is_valid_jacobian_ directly.
 * This causes Trust Region (and any solver that copies a group and
 * then calls applyJacobian on the copy) to throw.
 *
 * Fix: after the base-class operator= / copy-ctor, transfer ownership
 * based on is_valid_jacobian_ (a protected member) alone.
 * ========================================================================= */
class MBDynThyraGroup : public NOX::Thyra::Group
{
public:
     // Forward the "power user" constructor used by BuildSolver().
     MBDynThyraGroup(
          const NOX::Thyra::Vector& initial_guess,
          const Teuchos::RCP<const Thyra::ModelEvaluator<double>>& model,
          const Teuchos::RCP<Thyra::LinearOpBase<double>>& linear_op,
          const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>& lows_factory,
          const Teuchos::RCP<Thyra::PreconditionerBase<double>>& prec_op,
          const Teuchos::RCP<Thyra::PreconditionerFactoryBase<double>>& prec_factory)
          : NOX::Thyra::Group(initial_guess, model, linear_op,
                              lows_factory, prec_op, prec_factory)
     {}

     // Copy constructor (used by clone).
     MBDynThyraGroup(const MBDynThyraGroup& source, NOX::CopyType type)
          : NOX::Thyra::Group(source, type)
     {
          if (type == NOX::DeepCopy
              && Teuchos::nonnull(shared_jacobian_)
              && is_valid_jacobian_)
          {
               shared_jacobian_->getObject(this);
          }
     }

     NOX::Abstract::Group&
     operator=(const NOX::Abstract::Group& source) override
     {
          NOX::Thyra::Group::operator=(source);
          if (Teuchos::nonnull(shared_jacobian_) && is_valid_jacobian_)
               shared_jacobian_->getObject(this);
          return *this;
     }

     Teuchos::RCP<NOX::Abstract::Group>
     clone(NOX::CopyType type) const override
     {
          return Teuchos::rcp(new MBDynThyraGroup(*this, type));
     }

     // NOX::Thyra::Group does not implement computeGradient() (it
     // always throws).  The Epetra version computes grad = J^T * F.
     // Implement the same here so that line search methods like
     // More-Thuente that need the gradient work correctly.
     NOX::Abstract::Group::ReturnType
     computeGradient() override
     {
          if (isGradient())
               return NOX::Abstract::Group::Ok;

          if (!isF()) {
               silent_cerr("MBDynThyraGroup::computeGradient() - "
                           "residual (F) is not valid" << std::endl);
               return NOX::Abstract::Group::BadDependency;
          }

          if (!isJacobian()) {
               silent_cerr("MBDynThyraGroup::computeGradient() - "
                           "Jacobian is not valid" << std::endl);
               return NOX::Abstract::Group::BadDependency;
          }

          // gradient = J^T * F
          const auto& jac = shared_jacobian_->getObject(this);
          const auto& f_thyra = dynamic_cast<const NOX::Thyra::Vector&>(
               getF()).getThyraRCPVector();
          auto g_thyra = dynamic_cast<NOX::Thyra::Vector&>(
               *gradient_vec_).getThyraRCPVector();

          ::Thyra::apply(*jac, ::Thyra::TRANS, *f_thyra, g_thyra.ptr());

          is_valid_gradient_dir_ = true;
          return NOX::Abstract::Group::Ok;
     }
};

/* =========================================================================
 * NoxNonlinearSolver
 * ========================================================================= */
class NoxNonlinearSolver : public NonlinearSolver,
                           private NoxSolverParameters,
                           private NOX::Abstract::PrePostOperator
{
public:
     friend class ModelEvaluatorWrapper;
     friend class TpetraMatFreeJacOper;
     friend class MBDynJacobianOp;
     friend class MBDynLinearOpWithSolve;
     friend class MBDynPrecOp;
     friend class NoxResidualTest;
     friend class NoxSolutionTest;

     NoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                        const NoxSolverParameters& oParam);
     ~NoxNonlinearSolver();

     void Solve(const NonlinearProblem* pNLP,
                Solver* pS,
                const integer iMaxIter,
                const doublereal& Tol,
                integer& iIterCnt,
                doublereal& dErr,
                const doublereal& SolTol,
                doublereal& dSolErr) override;

     bool NoxMakeSolTest(const VectorHandler& XPrev,
                         const VectorHandler& XCurr,
                         const doublereal& dTol,
                         doublereal& dTest);

     bool NoxMakeResTest(const VectorHandler& oResVec,
                         const doublereal& dTol,
                         doublereal& dTest,
                         doublereal& dTestDiff);

private:
     struct CPUTimeGuard {
          explicit CPUTimeGuard(const NoxNonlinearSolver& s, CPUTimeType t)
               : oWatch(s, t) { oWatch.Tic(); }
          ~CPUTimeGuard() { oWatch.Toc(); }
          CPUStopWatch oWatch;
     };

     void Attach(Solver* pS, const NonlinearProblem* pNLP);
     void Residual(const VectorHandler* pSol, VectorHandler* pRes);
     void Jacobian();
     void BuildSolver(integer iMaxIter);
     void OutputIteration(integer iIterCnt, bool bJacobian) const;

     void runPreIterate(const NOX::Solver::Generic&)            override {}
     void runPostIterate(const NOX::Solver::Generic&)           override {}
     void runPreSolve(const NOX::Solver::Generic&)              override {}
     void runPostSolve(const NOX::Solver::Generic&)             override {}
     void runPreSolutionUpdate(const NOX::Abstract::Vector&,
                               const NOX::Solver::Generic&)     override {}
     void runPreLineSearch(const NOX::Solver::Generic& solver)  override;
     void runPostLineSearch(const NOX::Solver::Generic& solver) override;

     void ResetPrecondReuse() const {
          iPrecInnerIterCnt = iPrecInnerIterCntTot = 0;
     }
     void ForcePrecondRebuild() const {
          iPrecInnerIterCnt = iPrecInnerIterCntTot =
               std::numeric_limits<integer>::max();
     }

     Teuchos::RCP<const TpetraComm>      pComm;
     Teuchos::RCP<ModelEvaluatorWrapper> pModelEval;
     Teuchos::RCP<NOX::Solver::Generic>  pNonlinearSolver;
     Teuchos::RCP<NOX::Thyra::Vector>    pSolutionView;
     Teuchos::ParameterList              oSolverParam;

     const NonlinearProblem* pNonlinearProblem;
     Solver*                 pSolver;
     SolutionManager*        pSolutionManager;

     mutable MyVectorHandler DeltaX, XPrev, TmpRes;

     NoxResidualTest oResTest;
     NoxSolutionTest oSolTest;

     mutable bool bUseTranspose;
     bool         bUpdateJacobian;
     bool         bInDerivativeSolver;
     bool         bInLineSearch;

     mutable integer iPrecInnerIterCnt;
     mutable integer iPrecInnerIterCntTot;
     mutable integer iInnerIterCntTot;
};

/* =========================================================================
 * TpetraMatFreeJacOper::applyImpl
 * ========================================================================= */
void TpetraMatFreeJacOper::applyImpl(
     const Thyra::EOpTransp                          M_trans,
     const Thyra::MultiVectorBase<SC>&               X_in,
     const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
     const SC                                        alpha,
     const SC                                        beta) const
{
     if (M_trans != Thyra::NOTRANS) {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
     // NOX always calls with alpha=1, beta=0
     ASSERT(alpha == SC(1.) && beta == SC(0.));

     NoxNonlinearSolver::CPUTimeGuard oCPUTimeJac(
          oNoxSolver, NoxNonlinearSolver::CPU_JACOBIAN);

     ASSERT(oNoxSolver.Size > 0);

     auto xTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraMultiVector(Teuchos::rcpFromRef(X_in));
     auto yTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getTpetraMultiVector(Teuchos::rcpFromRef(*Y_out));

     ASSERT(static_cast<integer>(xTpetra->getLocalLength()) == oNoxSolver.Size);
     ASSERT(static_cast<integer>(yTpetra->getLocalLength()) == oNoxSolver.Size);

     const MyVectorHandler X(
          oNoxSolver.Size,
          const_cast<doublereal*>(
               xTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data()));
     MyVectorHandler Y(
          oNoxSolver.Size,
          yTpetra->getLocalViewHost(Tpetra::Access::ReadWrite).data());

     // Matrix-free product: Y = J * X
     oNoxSolver.pNonlinearProblem->Jacobian(&Y, &X);

#ifdef DEBUG_JACOBIAN
     if (!pA) {
          SAFENEWWITHCONSTRUCTOR(pA,
                                  SpGradientSparseMatrixHandler,
                                  SpGradientSparseMatrixHandler(
                                       oNoxSolver.Size, oNoxSolver.Size));
          oNoxSolver.pNonlinearProblem->Jacobian(pA);
          pA->PacMat();
          AX.Resize(oNoxSolver.Size);
     }
     pA->MatVecMul(AX, X);
     const doublereal dTolRel = sqrt(std::numeric_limits<doublereal>::epsilon());
     const doublereal dTolAbs = sqrt(std::numeric_limits<doublereal>::epsilon());
     const doublereal dNormAX = AX.Norm();
     for (integer i = 1; i <= X.iGetSize(); ++i) {
          if (std::fabs(AX(i) - Y(i)) > dTolRel + dTolAbs * dNormAX) {
               DEBUGCERR("MatFree Jacobian check failed: AX("
                         << i << ")=" << AX(i)
                         << " Y(" << i << ")=" << Y(i) << "\n");
               ASSERT(0);
          }
     }
#endif
}

/* =========================================================================
 * MBDynPrecOp::applyImpl  out-of-line implementation
 *
 * Applies the preconditioner M^{-1} * x via pSolutionManager->Solve().
 * This mirrors NoxNonlinearSolver::ApplyInverse from the Epetra version.
 * ========================================================================= */
void MBDynPrecOp::applyImpl(
     const Thyra::EOpTransp                          M_trans,
     const Thyra::MultiVectorBase<SC>&               X_in,
     const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
     const SC                                        alpha,
     const SC                                        beta) const
{
     ASSERT(M_trans == Thyra::NOTRANS);
     ASSERT(alpha == SC(1.) && beta == SC(0.));
     ASSERT(oNoxSolver.pSolutionManager != nullptr);

     ++oNoxSolver.iPrecInnerIterCnt;
     ++oNoxSolver.iInnerIterCntTot;

     NoxNonlinearSolver::CPUTimeGuard oCPULinearSolver(
          oNoxSolver, NoxNonlinearSolver::CPU_LINEAR_SOLVER);

     VectorHandler* const pResVec = oNoxSolver.pSolutionManager->pResHdl();
     const VectorHandler* const pSolVec = oNoxSolver.pSolutionManager->pSolHdl();
     const integer n = oNoxSolver.Size;

     ASSERT(n == pResVec->iGetSize());
     ASSERT(n == pSolVec->iGetSize());

     auto bTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraMultiVector(Teuchos::rcpFromRef(X_in));
     auto xTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getTpetraMultiVector(Teuchos::rcpFromRef(*Y_out));

     const SC* bData =
          bTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data();
     SC* xData =
          xTpetra->getLocalViewHost(Tpetra::Access::ReadWrite).data();

     std::copy(bData, bData + n, pResVec->pdGetVec());

     oNoxSolver.pSolutionManager->Solve();

     std::copy(pSolVec->pdGetVec(), pSolVec->pdGetVec() + n, xData);
}

/* =========================================================================
 * MBDynJacobianOp  out-of-line implementation
 * (must follow full definition of NoxNonlinearSolver)
 * ========================================================================= */
void MBDynJacobianOp::applyImpl(
     const Thyra::EOpTransp                          M_trans,
     const Thyra::MultiVectorBase<SC>&               X_in,
     const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
     const SC                                        alpha,
     const SC                                        beta) const
{
     ASSERT(alpha == SC(1.) && beta == SC(0.));
     ASSERT(oNoxSolver.pSolutionManager != nullptr);

     const MatrixHandler* const pJacMat =
          oNoxSolver.pSolutionManager->pMatHdl();
     ASSERT(pJacMat != nullptr);

     auto xTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraMultiVector(Teuchos::rcpFromRef(X_in));
     auto yTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getTpetraMultiVector(Teuchos::rcpFromRef(*Y_out));

     const integer n = oNoxSolver.Size;
     const MyVectorHandler XVec(
          n,
          const_cast<SC*>(
               xTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data()));
     MyVectorHandler YVec(
          n,
          yTpetra->getLocalViewHost(Tpetra::Access::ReadWrite).data());

     if (M_trans == Thyra::TRANS) {
          pJacMat->MatTVecMul(YVec, XVec);
     } else {
          pJacMat->MatVecMul(YVec, XVec);
     }
}

/* =========================================================================
 * MBDynLinearOpWithSolve  out-of-line implementation
 * (must follow full definition of NoxNonlinearSolver)
 * ========================================================================= */
void MBDynLinearOpWithSolve::applyImpl(
     const Thyra::EOpTransp                          M_trans,
     const Thyra::MultiVectorBase<SC>&               X_in,
     const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
     const SC                                        alpha,
     const SC                                        beta) const
{
     ASSERT(alpha == SC(1.) && beta == SC(0.));
     ASSERT(oNoxSolver.pSolutionManager != nullptr);

     const MatrixHandler* const pJacMat =
          oNoxSolver.pSolutionManager->pMatHdl();
     ASSERT(pJacMat != nullptr);

     auto xTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraMultiVector(Teuchos::rcpFromRef(X_in));
     auto yTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getTpetraMultiVector(Teuchos::rcpFromRef(*Y_out));

     const integer n = oNoxSolver.Size;
     const MyVectorHandler XVec(
          n,
          const_cast<SC*>(
               xTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data()));
     MyVectorHandler YVec(
          n,
          yTpetra->getLocalViewHost(Tpetra::Access::ReadWrite).data());

     if (M_trans == Thyra::TRANS) {
          pJacMat->MatTVecMul(YVec, XVec);
     } else {
          pJacMat->MatVecMul(YVec, XVec);
     }
}

Thyra::SolveStatus<TpetraSC> MBDynLinearOpWithSolve::solveImpl(
     const Thyra::EOpTransp                          M_trans,
     const Thyra::MultiVectorBase<SC>&               B,
     const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& X,
     Teuchos::Ptr<const Thyra::SolveCriteria<SC>>    /* solveCriteria */) const
{
     ASSERT(M_trans == Thyra::NOTRANS);
     ASSERT(oNoxSolver.pSolutionManager != nullptr);

     ++oNoxSolver.iPrecInnerIterCnt;
     ++oNoxSolver.iInnerIterCntTot;

     NoxNonlinearSolver::CPUTimeGuard oCPULinearSolver(
          oNoxSolver, NoxNonlinearSolver::CPU_LINEAR_SOLVER);

     VectorHandler* const pResVec = oNoxSolver.pSolutionManager->pResHdl();
     const VectorHandler* const pSolVec = oNoxSolver.pSolutionManager->pSolHdl();
     const integer n = oNoxSolver.Size;

     ASSERT(n == pResVec->iGetSize());
     ASSERT(n == pSolVec->iGetSize());

     auto bTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraMultiVector(Teuchos::rcpFromRef(B));
     auto xTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getTpetraMultiVector(Teuchos::rcpFromRef(*X));

     const SC* bData =
          bTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data();
     SC* xData =
          xTpetra->getLocalViewHost(Tpetra::Access::ReadWrite).data();

     // Copy RHS into MBDyn's residual vector (mirrors old ApplyInverse)
     std::copy(bData, bData + n, pResVec->pdGetVec());

     oNoxSolver.pSolutionManager->Solve();

     // Copy MBDyn's solution back into the Thyra output vector
     std::copy(pSolVec->pdGetVec(), pSolVec->pdGetVec() + n, xData);

     Thyra::SolveStatus<SC> solveStatus;
     solveStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
     return solveStatus;
}

/* =========================================================================
 * ModelEvaluatorWrapper  implementation
 * ========================================================================= */
ModelEvaluatorWrapper::ModelEvaluatorWrapper(
     NoxNonlinearSolver& solver_a,
     integer iSize,
     const Teuchos::RCP<const TpetraComm>& pComm_a)
     : oNoxSolver(solver_a), pComm(pComm_a)
{
     Rebuild(iSize);
}

void ModelEvaluatorWrapper::Rebuild(integer iSize)
{
     pMap = Teuchos::rcp(new TpetraMap(
          static_cast<TpetraGO>(iSize), TpetraGO(0), pComm));
     pSpace =
          Thyra::createVectorSpace<TpetraSC, TpetraLO, TpetraGO, TpetraNode>(pMap);

     pJacobianOp = Teuchos::null; // invalidate on resize
     pPreconditioner = Teuchos::null;

     // For the explicit-matrix path use MBDynLOWSFactory so the linear solve
     // is routed through pSolutionManager->Solve() regardless of which
     // concrete MatrixHandler / SolutionManager MBDyn is using.
     // For the matrix-free path pLOWSFactory is set by BuildSolver() via
     // Stratimikos after this call returns.
     if (!(oNoxSolver.uFlags & NoxSolverParameters::JACOBIAN_NEWTON_KRYLOV)
         || (oNoxSolver.uFlags & NoxSolverParameters::USE_PRECOND_AS_SOLVER)) {
          // For the explicit-matrix path, or JFNK with "use preconditioner
          // as solver", the linear solve is routed through
          // pSolutionManager->Solve().
          pLOWSFactory = Teuchos::rcp(new MBDynLOWSFactory(oNoxSolver, pSpace));
     }

     oNominalValues = this->createInArgs();
     auto x0 = Thyra::createMember(pSpace);
     Thyra::assign(x0.ptr(), TpetraSC(0.));
     oNominalValues.set_x(x0);
}

Thyra::ModelEvaluatorBase::InArgs<TpetraSC>
ModelEvaluatorWrapper::createInArgs() const
{
     Thyra::ModelEvaluatorBase::InArgsSetup<SC> inArgs;
     inArgs.setModelEvalDescription(this->description());
     inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x);
     return inArgs;
}

Thyra::ModelEvaluatorBase::OutArgs<TpetraSC>
ModelEvaluatorWrapper::createOutArgsImpl() const
{
     Thyra::ModelEvaluatorBase::OutArgsSetup<SC> outArgs;
     outArgs.setModelEvalDescription(this->description());
     outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
     outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
     outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec, true);
     return outArgs;
}

Teuchos::RCP<Thyra::LinearOpBase<TpetraSC>>
ModelEvaluatorWrapper::create_W_op() const
{
     // Lazily construct and cache.
     // pSolutionManager must already be valid (Attach() precedes BuildSolver()).
     if (pJacobianOp.is_null()) {
          if (oNoxSolver.uFlags & NoxSolverParameters::JACOBIAN_NEWTON_KRYLOV) {
               // Matrix-free: recomputes J*v on every apply call.
               pJacobianOp = Teuchos::rcp(
                    new TpetraMatFreeJacOper(oNoxSolver, pSpace));
          } else {
               // Explicit matrix: wrap MBDyn's MatrixHandler generically.
               // MBDynJacobianOp calls pMatHdl()->MatVecMul() which works with
               // *any* matrix type that MBDyn's SolutionManager exposes.
               // This replaces the old TpetraLinearOp / dynamic_cast path that
               // assumed pMatHdl() returns a TpetraSparseMatrixHandler.
               ASSERT(oNoxSolver.pSolutionManager != nullptr);
               pJacobianOp = Teuchos::rcp(
                    new MBDynJacobianOp(oNoxSolver, pSpace));
          }
     }
     return pJacobianOp;
}

Teuchos::RCP<Thyra::PreconditionerBase<TpetraSC>>
ModelEvaluatorWrapper::create_W_prec() const
{
     if (pPreconditioner.is_null()) {
          ASSERT(oNoxSolver.pSolutionManager != nullptr);
          auto pPrecOp = Teuchos::rcp(
               new MBDynPrecOp(oNoxSolver, pSpace));
          Teuchos::RCP<Thyra::LinearOpBase<SC>> pPrecOpBase = pPrecOp;
          // Wrap as "unspecified" (Belos will use it as right prec).
          pPreconditioner = Teuchos::rcp(
               new Thyra::DefaultPreconditioner<SC>(pPrecOpBase));
     }
     return pPreconditioner;
}

void ModelEvaluatorWrapper::evalModelImpl(
     const Thyra::ModelEvaluatorBase::InArgs<SC>&  inArgs,
     const Thyra::ModelEvaluatorBase::OutArgs<SC>& outArgs) const
{
     /* --- Unpack x (zero-copy) ------------------------------------------ */
     auto xThyra = inArgs.get_x();
     ASSERT(!xThyra.is_null());

     auto xTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraVector(xThyra);

     const MyVectorHandler oSol(
          oNoxSolver.Size,
          const_cast<doublereal*>(
               xTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data()));

     /* --- Residual branch  (mirrors computeF) --------------------------- */
     if (!outArgs.get_f().is_null()) {
          auto fThyra  = outArgs.get_f();
          auto fTpetra =
               Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
                    ::getTpetraVector(fThyra);

          MyVectorHandler oRes(
               oNoxSolver.Size,
               fTpetra->getLocalViewHost(Tpetra::Access::ReadWrite).data());

          oNoxSolver.bUpdateJacobian = true;

          if (oNoxSolver.bInLineSearch && oNoxSolver.pNonlinearSolver) {
               const auto& oLS =
                    dynamic_cast<const NOX::Solver::LineSearchBased&>(
                         *oNoxSolver.pNonlinearSolver);
               oNoxSolver.SetNonlinearSolverHint(
                    NonlinearSolver::LINESEARCH_LAMBDA_CURR,
                    oLS.getStepSize());
               DEBUGCERR("line search iteration "
                    << oNoxSolver.GetNonlinearSolverHint(
                         NonlinearSolver::LINESEARCH_ITERATION_CURR)
                    << ": lambda=" << oLS.getStepSize() << "\n");
          }

          oNoxSolver.Residual(&oSol, &oRes);
          // std::cerr << "chiamato Residual" << std::endl;
          //           std::cerr << "==============" << std::endl;
          //           oNoxSolver.pSolver->PrintSolution(
          //                oNoxSolver.DeltaX,
          //                oNoxSolver.pNonlinearSolver->getNumIterations());
          //           std::cerr << "--------------" << std::endl;
          //           std::cerr << oSol << std::endl;
          //           std::cerr << ";;;;;;;;;;;;;;" << std::endl;
          //           oNoxSolver.pSolver->PrintResidual(
          //                oRes,
          //                oNoxSolver.pNonlinearSolver->getNumIterations());
          //           std::cerr << "**************" << std::endl;
          //           std::cerr << oRes << std::endl;
          //           std::cerr << "++++++++++++++" << std::endl;

          if (oNoxSolver.pSolver && oNoxSolver.pNonlinearSolver) {
               if (oNoxSolver.outputSol()) {
                    oNoxSolver.pSolver->PrintSolution(
                         oNoxSolver.DeltaX,
                         oNoxSolver.pNonlinearSolver->getNumIterations());
               }
               if (oNoxSolver.outputRes()) {
                    oNoxSolver.pSolver->PrintResidual(
                         oRes,
                         oNoxSolver.pNonlinearSolver->getNumIterations());
               }
          }

          oRes *= -1.; // NOX sign convention

          if (oNoxSolver.bInLineSearch) {
               const integer iIterCurr =
                    oNoxSolver.GetNonlinearSolverHint(
                         NonlinearSolver::LINESEARCH_ITERATION_CURR);
               oNoxSolver.SetNonlinearSolverHint(
                    NonlinearSolver::LINESEARCH_ITERATION_CURR,
                    iIterCurr + 1);
          }
     }

     /* --- Jacobian branch  (mirrors computeJacobian) -------------------- */
     // For the explicit-matrix path, assembling below into
     // pSolutionManager->pMatHdl() automatically updates the operator NOX
     // uses for the linear solve — no copy required.
     // For the matrix-free (JFNK) path, only the MBDyn state needs to be
     // current; no explicit matrix assembly is performed (the matrix-free
     // operator recomputes J*v on every apply).
     if (!outArgs.get_W_op().is_null()) {
          if (oNoxSolver.bUpdateJacobian) {
               // MBDyn convention: Residual() must precede Jacobian().
               // If the residual branch above already ran, the solution is
               // current; otherwise call explicitly.
               if (outArgs.get_f().is_null()) {
                    oNoxSolver.Residual(&oSol, &oNoxSolver.TmpRes);
               }
               if (!(oNoxSolver.uFlags
                     & NoxSolverParameters::JACOBIAN_NEWTON_KRYLOV)) {
                    // Explicit-matrix path: assemble the Jacobian matrix.
                    if (oNoxSolver.pSolutionManager) {
                         oNoxSolver.pSolutionManager->MatrReset();
                    }
                    oNoxSolver.Jacobian();
                    // Clear the flag: the assembled matrix is now current.
                    oNoxSolver.bUpdateJacobian = false;
               }
               // For JFNK: do NOT clear bUpdateJacobian here.
               // The flag is consumed by the W_prec branch (separate
               // evalModel call) which needs to know it should assemble
               // the preconditioner matrix.
          }
     }

     /* --- Preconditioner branch  (mirrors computePreconditioner) --------- */
     // For JFNK, assemble the explicit Jacobian matrix for use as
     // preconditioner.  The MBDynPrecOp operator reads from
     // pSolutionManager->pMatHdl() and applies J^{-1} via Solve().
     //
     // In the Epetra version, the NOX::Epetra::LinearSystem called
     // getPreconditionerPolicy() before each linear solve to decide
     // whether to recompute the preconditioner.  The policy checked:
     //   iPrecInnerIterCnt >= iInnerIterBeforeAssembly  (GMRES iters in prev step)
     //   iPrecInnerIterCntTot >= iIterationsBeforeAssembly (total across Newton steps)
     // If either was true → RECOMPUTE, else → REUSE.
     // Then recomputePreconditioner → computeJacobian (guarded by bUpdateJacobian).
     //
     // In the Thyra port, NOX::Thyra::Group calls evalModel(W_prec)
     // every Newton iteration (from updateLOWS).  We replicate the
     // Epetra reuse policy here: only assemble when ForcePrecondRebuild
     // was called or the inner-iteration count exceeds the threshold.
     if (!outArgs.get_W_prec().is_null()) {
          const bool bPrecRecompute =
               oNoxSolver.iPrecInnerIterCnt >= oNoxSolver.iInnerIterBeforeAssembly
               || oNoxSolver.iPrecInnerIterCntTot >= oNoxSolver.iIterationsBeforeAssembly;

          if (bPrecRecompute && oNoxSolver.bUpdateJacobian) {
               // MBDyn convention: Residual() must precede Jacobian().
               if (outArgs.get_f().is_null()) {
                    oNoxSolver.Residual(&oSol, &oNoxSolver.TmpRes);
               }
               if (oNoxSolver.pSolutionManager) {
                    oNoxSolver.pSolutionManager->MatrReset();
               }
               oNoxSolver.Jacobian();
               oNoxSolver.bUpdateJacobian = false;
               oNoxSolver.ResetPrecondReuse();
          } else if (!bPrecRecompute) {
               ++oNoxSolver.iPrecInnerIterCntTot;
          }
          // Reset per-Newton-step GMRES counter so it only counts the
          // next linear solve's preconditioner applications.  Mirrors
          // iPrecInnerIterCnt = 0 at the start of applyJacobianInverse
          // in the Epetra version.  (ResetPrecondReuse already cleared
          // it in the recompute path, so this only matters for reuse.)
          oNoxSolver.iPrecInnerIterCnt = 0;
     }
}

/* =========================================================================
 * NoxResidualTest  implementation
 * ========================================================================= */
NOX::StatusTest::StatusType
NoxResidualTest::checkStatus(const NOX::Solver::Generic& problem,
                              NOX::StatusTest::CheckType  checkType)
{
     if (checkType == NOX::StatusTest::None) {
          eStatus = NOX::StatusTest::Unevaluated;
          return eStatus;
     }
     const NOX::Abstract::Group& grp = problem.getSolutionGroup();
     if (!grp.isF()) {
          eStatus = NOX::StatusTest::Unevaluated;
          return eStatus;
     }

     const auto& FT =
          dynamic_cast<const NOX::Thyra::Vector&>(grp.getF());
     auto fTpetra =
          Thyra::TpetraOperatorVectorExtraction<TpetraSC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraVector(
                    Teuchos::rcpFromRef(FT.getThyraVector()));

     const MyVectorHandler oResVec(
          static_cast<integer>(fTpetra->getLocalLength()),
          const_cast<doublereal*>(
               fTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data()));

     const doublereal dFirstResFact =
          problem.getNumIterations() == 0 ? 1e-2 : 1.;

     eStatus = oNoxSolver.NoxMakeResTest(
          oResVec, dFirstResFact * dTolRes, dErrRes, dErrResDiff)
          ? NOX::StatusTest::Converged
          : NOX::StatusTest::Unconverged;
     return eStatus;
}

std::ostream& NoxResidualTest::print(std::ostream& stream, int indent) const
{
     for (int j = 0; j < indent; ++j) stream << ' ';
     stream << eStatus
            << "F-Norm = " << NOX::Utils::sciformat(dErrRes, 3)
            << " < "       << NOX::Utils::sciformat(dTolRes, 3) << "\n";
     return stream;
}

/* =========================================================================
 * NoxSolutionTest  implementation
 * ========================================================================= */
NOX::StatusTest::StatusType
NoxSolutionTest::checkStatus(const NOX::Solver::Generic& problem,
                              NOX::StatusTest::CheckType  checkType)
{
     if (checkType == NOX::StatusTest::None) {
          eStatus = NOX::StatusTest::Unevaluated;
          return eStatus;
     }
     if (problem.getNumIterations() == 0) {
          eStatus = NOX::StatusTest::Unconverged;
          return eStatus;
     }
     if (!problem.getSolutionGroup().isF()) {
          eStatus = NOX::StatusTest::Unconverged;
          return eStatus;
     }

     const integer n = oNoxSolver.Size;
     auto extractRaw = [n](const NOX::Abstract::Vector& VA) {
          const auto& VT = dynamic_cast<const NOX::Thyra::Vector&>(VA);
          auto tp =
               Thyra::TpetraOperatorVectorExtraction<
                    TpetraSC, TpetraLO, TpetraGO, TpetraNode>
                         ::getConstTpetraVector(
                              Teuchos::rcpFromRef(VT.getThyraVector()));
          return MyVectorHandler(
               n,
               const_cast<doublereal*>(
                    tp->getLocalViewHost(Tpetra::Access::ReadOnly).data()));
     };

     MyVectorHandler XPrevVec =
          extractRaw(problem.getPreviousSolutionGroup().getX());
     MyVectorHandler XCurrVec =
          extractRaw(problem.getSolutionGroup().getX());

     eStatus = oNoxSolver.NoxMakeSolTest(XPrevVec, XCurrVec, dTolSol, dErrSol)
          ? NOX::StatusTest::Converged
          : NOX::StatusTest::Unconverged;
     return eStatus;
}

std::ostream& NoxSolutionTest::print(std::ostream& stream, int indent) const
{
     for (int j = 0; j < indent; ++j) stream << ' ';
     stream << eStatus
            << "X-Norm = " << NOX::Utils::sciformat(dErrSol, 3)
            << " < "       << NOX::Utils::sciformat(dTolSol, 3) << "\n";
     return stream;
}

/* =========================================================================
 * NoxNonlinearSolver  implementation
 * ========================================================================= */
NoxNonlinearSolver::NoxNonlinearSolver(
     const NonlinearSolverTestOptions& oSolverOpt,
     const NoxSolverParameters& oParam)
     : NonlinearSolver(oSolverOpt),
       NoxSolverParameters(oParam),
       pNonlinearProblem(nullptr),
       pSolver(nullptr),
       pSolutionManager(nullptr),
       oResTest(*this),
       oSolTest(*this),
       bUseTranspose(false),
       bUpdateJacobian(true),
       bInDerivativeSolver(true),
       bInLineSearch(false),
       iInnerIterCntTot(0)
{
#ifdef USE_MPI
     pComm = tpetraMpiComm(MBDynComm);
#else
     pComm = tpetraSerialComm();
#endif
     ForcePrecondRebuild();
}

NoxNonlinearSolver::~NoxNonlinearSolver()
{
     silent_cerr("total inner iterations: " << iInnerIterCntTot << "\n");

     // Release Trilinos objects in dependency order BEFORE the
     // SolutionManager (owned by Solver) is destroyed.  The NOX
     // solver, model evaluator and solution view hold RCPs to
     // Tpetra operators/vectors that reference the SolutionManager's
     // CrsMatrix and MultiVectors via raw pointers.
     pNonlinearSolver.reset();
     pSolutionView.reset();
     pModelEval.reset();
     pComm.reset();
}

void
NoxNonlinearSolver::Solve(const NonlinearProblem* pNLP,
                           Solver* pS,
                           const integer iMaxIter,
                           const doublereal& dTolRes,
                           integer& iIterCnt,
                           doublereal& dResErr,
                           const doublereal& dTolSol,
                           doublereal& dSolErr)
{
     DEBUGCERR("Solve()\n");

     bInLineSearch = false;
     SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, 0);
     SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.);

     oResTest.SetTolerance(dTolRes);
     oSolTest.SetTolerance(dTolSol);

     Attach(pS, pNLP);

     if (!pNonlinearSolver.get()) {
          BuildSolver(iMaxIter);
     }

     pNonlinearSolver->reset(*pSolutionView);

     if (!bKeepJacAcrossSteps || bInDerivativeSolver) {
          ForcePrecondRebuild();
     }

     iIterCnt = 0;

     for (;;) {
          oResTest.Reset();
          oSolTest.Reset();

          const integer TotJacPrev = TotJac;

          DEBUGCERR("Nonlinear solver step(" << iIterCnt << ")\n");

          NOX::StatusTest::StatusType solvStatus = pNonlinearSolver->step();

          if (oResTest.getStatus() == NOX::StatusTest::Unevaluated) {
               oResTest.checkStatus(*pNonlinearSolver,
                                    NOX::StatusTest::CheckType::Complete);
          }

          ++iIterCnt;

          OutputIteration(iIterCnt, TotJac > TotJacPrev);
          pSolver->CheckTimeStepLimit(oResTest.dGetTest(),
                                       oResTest.dGetTestDiff());

          dResErr = oResTest.dGetTest();
          dSolErr = (oSolTest.getStatus() != NOX::StatusTest::Unevaluated)
               ? oSolTest.dGetTest() : -1.;

          if (solvStatus == NOX::StatusTest::Converged) break;

          if (iIterCnt > iMaxIter || solvStatus == NOX::StatusTest::Failed) {
               throw NoConvergence(MBDYN_EXCEPT_ARGS);
          }

          if (mbdyn_stop_at_end_of_iteration()) {
               throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
          }
     }
}

bool NoxNonlinearSolver::NoxMakeSolTest(const VectorHandler& XPrev_a,
                                         const VectorHandler& XCurr,
                                         const doublereal& dTol,
                                         doublereal& dTest)
{
     DeltaX.ScalarAddMul(XCurr, XPrev_a, -1.);
     return NonlinearSolver::MakeSolTest(pSolver, DeltaX, dTol, dTest);
}

bool NoxNonlinearSolver::NoxMakeResTest(const VectorHandler& oResVec,
                                         const doublereal& dTol,
                                         doublereal& dTest,
                                         doublereal& dTestDiff)
{
     return NonlinearSolver::MakeResTest(pSolver, pNonlinearProblem,
                                          oResVec, dTol, dTest, dTestDiff);
}

void NoxNonlinearSolver::BuildSolver(const integer iMaxIter_a)
{
     // Attach() must have been called first so pSolutionManager is valid
     // when create_W_op() resolves the operator.
     if (!pModelEval) {
          pModelEval = Teuchos::rcp(
               new ModelEvaluatorWrapper(*this, Size, pComm));
     } else {
          pModelEval->Rebuild(Size);
     }

     auto x0 = Thyra::createMember(pModelEval->GetSpace());
     Thyra::assign(x0.ptr(), TpetraSC(0.));
     pSolutionView = Teuchos::rcp(new NOX::Thyra::Vector(x0));

     /* ---- Linear solve factory ----------------------------------------- */
     // For the explicit-matrix path (any SolutionManager):
     //   MBDynLOWSFactory was already created in ModelEvaluatorWrapper::Rebuild()
     //   and delegates to pSolutionManager->Solve() — works with any linear solver
     //   that MBDyn wraps (Umfpack, SuperLU, PARDISO, KLU, Tpetra-based, etc.).
     //
     // For the matrix-free path:
     //   Stratimikos/Belos provides the Krylov solver; no preconditioner is set
     //   (same as original Epetra path with USE_PRECOND_AS_SOLVER not set).
     Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<TpetraSC>> pLOWSFactory;

     if (uFlags & JACOBIAN_NEWTON_KRYLOV) {
          if (uFlags & USE_PRECOND_AS_SOLVER) {
               // When "use preconditioner as solver" is enabled, bypass
               // the Krylov iteration entirely and route the linear solve
               // through pSolutionManager->Solve() (the direct solver /
               // preconditioner).  This mirrors the Epetra code path that
               // called ApplyInverse() directly in applyJacobianInverse().
               pLOWSFactory = pModelEval->GetLOWSFactory();
               ASSERT(!pLOWSFactory.is_null());
          } else {
          /* ---- Stratimikos / Belos (matrix-free Krylov) ------------------ */
          Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
          Teuchos::RCP<Teuchos::ParameterList> pLSParams =
               Teuchos::rcp(new Teuchos::ParameterList());

          pLSParams->set("Linear Solver Type", "Belos");
          pLSParams->set("Preconditioner Type", "None");

          std::string sBelosSolverType = "Block GMRES";
          if      (uFlags & LINEAR_SOLVER_PSEUDO_BLOCK_GMRES)         sBelosSolverType = "Pseudo Block GMRES";
          else if (uFlags & LINEAR_SOLVER_BLOCK_CG)                   sBelosSolverType = "Block CG";
          else if (uFlags & LINEAR_SOLVER_PSEUDO_BLOCK_CG)            sBelosSolverType = "Pseudo Block CG";
          else if (uFlags & LINEAR_SOLVER_BLOCK_STOCHASTIC_CG)        sBelosSolverType = "Pseudo Block Stochastic CG";
          else if (uFlags & LINEAR_SOLVER_GCRODR)                     sBelosSolverType = "GCRODR";
          else if (uFlags & LINEAR_SOLVER_RCG)                        sBelosSolverType = "RCG";
          else if (uFlags & LINEAR_SOLVER_MINRES)                     sBelosSolverType = "MINRES";
          else if (uFlags & LINEAR_SOLVER_TFQMR)                      sBelosSolverType = "TFQMR";
          else if (uFlags & LINEAR_SOLVER_BICGSTAB)                   sBelosSolverType = "BiCGStab";
          else if (uFlags & LINEAR_SOLVER_FIXED_POINT)                sBelosSolverType = "Fixed Point";
          else if (uFlags & LINEAR_SOLVER_TPETRA_GMRES)               sBelosSolverType = "TPETRA GMRES";
          else if (uFlags & LINEAR_SOLVER_TPETRA_GMRES_PIPELINE)      sBelosSolverType = "TPETRA GMRES PIPELINE";
          else if (uFlags & LINEAR_SOLVER_TPETRA_GMRES_SINGLE_REDUCE) sBelosSolverType = "TPETRA GMRES SINGLE REDUCE";
          else if (uFlags & LINEAR_SOLVER_TPETRA_GMRES_SSTEP)         sBelosSolverType = "TPETRA GMRES S-STEP";

          auto& belosParams =
               pLSParams->sublist("Linear Solver Types").sublist("Belos");
          belosParams.set("Solver Type", sBelosSolverType);
          auto& bSolverParams =
               belosParams.sublist("Solver Types").sublist(sBelosSolverType);

          // ---- Obtain the authoritative valid-parameter list ---------------
          // eGetBelosParamType() already creates a throw-away solver for
          // this solver-type name on first call and caches its
          // getValidParameters() list; reuse that cache here rather than
          // duplicating the factory construction (which would also require
          // MV/OP/BelosSolver aliases that are not in scope in this TU).
          // UNKNOWN means "parameter absent from this solver's valid list".

          // ---- Named MBDyn parameters (always present in every solver) -----
          // "Maximum Iterations" and "Convergence Tolerance" appear in every
          // Belos solver's valid-parameter list.
          bSolverParams.set("Maximum Iterations",    iMaxIterLinSol);
          bSolverParams.set("Convergence Tolerance", dTolLinSol);

          // "Num Blocks" (Krylov restart depth) is only valid for GMRES-family
          // and the two recycling solvers.  Query the valid list so that the
          // code remains correct if Belos ever adds or removes the parameter
          // for a given solver, instead of relying on a hand-maintained mask.
          if (eGetBelosParamType(sBelosSolverType, "Num Blocks") != BelosParamType::UNKNOWN) {
               const integer iKrylovRestart =
                    std::min(Size, std::min(iMaxIterLinSol, iKrylovSubSpaceSize));
               bSolverParams.set("Num Blocks", iKrylovRestart);
          }

          bSolverParams.set("Output Frequency",
                            (uFlags & PRINT_CONVERGENCE_INFO) ? 1 : 0);
          bSolverParams.set("Output Style", 1);
          bSolverParams.set("Verbosity",
                            (uFlags & PRINT_CONVERGENCE_INFO)
                                 ? (Belos::Errors | Belos::Warnings |
                                    Belos::IterationDetails |
                                    Belos::StatusTestDetails)
                                 : Belos::Errors);

          // ---- Generic Belos parameter overrides from the input file -------
          // The "belos parameters" keyword block in the .mbd file is parsed
          // into oBelosParams as typed variants.  We apply each entry here
          // after confirming it is in the valid list.  Parameters set above
          // (Maximum Iterations, Convergence Tolerance, Num Blocks, Verbosity,
          // Output Frequency, Output Style) are intentionally overwritten only
          // by the named MBDyn keywords and are excluded from this loop via the
          // MBDYN_RESERVED_BELOS_PARAMS guard so that the solver's behaviour
          // remains consistent with what the MBDyn input says.
          static const std::set<std::string> MBDYN_RESERVED_BELOS_PARAMS = {
               "Maximum Iterations",
               "Convergence Tolerance",
               "Num Blocks",
               "Output Frequency",
               "Output Style",
               "Verbosity",
          };

          for (const auto& [sName, oValue] : oBelosParams) {
               if (MBDYN_RESERVED_BELOS_PARAMS.count(sName)) {
                    silent_cerr("Belos parameter \"" << sName
                                << "\" is controlled by a named MBDyn keyword "
                                << "and cannot be overridden in the "
                                << "\"belos parameters\" block; ignoring.\n");
                    continue;
               }
               if (eGetBelosParamType(sBelosSolverType, sName) == BelosParamType::UNKNOWN) {
                    silent_cerr("Belos parameter \"" << sName
                                << "\" is not valid for solver \""
                                << sBelosSolverType
                                << "\"; ignoring.\n");
                    continue;
               }
               // Type-safe dispatch using std::visit over the variant.
               std::visit([&bSolverParams, &sName](const auto& v) {
                    bSolverParams.set(sName, v);
               }, oValue);
          }

          linearSolverBuilder.setParameterList(pLSParams);
          pLOWSFactory = linearSolverBuilder.createLinearSolveStrategy("");
          // Store in pModelEval so NOX::Thyra::Group can retrieve it
          pModelEval->pLOWSFactory = pLOWSFactory;
          } // !USE_PRECOND_AS_SOLVER
     } else {
          // Explicit-matrix path: use the MBDynLOWSFactory already set in Rebuild()
          pLOWSFactory = pModelEval->GetLOWSFactory();
          ASSERT(!pLOWSFactory.is_null());
     }

     /* ---- NOX parameters ----------------------------------------------- */
     auto& oPrintParam    = oSolverParam.sublist("Printing");
     auto& oDirectionParam = oSolverParam.sublist("Direction");
     auto& oNewtonParam    = oDirectionParam.sublist("Newton");

     oNewtonParam.set("Forcing Term Minimum Tolerance", dForcingTermMinTol);
     oNewtonParam.set("Forcing Term Maximum Tolerance", dForcingTermMaxTol);
     oNewtonParam.set("Forcing Term Alpha",              dForcingTermAlpha);
     oNewtonParam.set("Forcing Term Gamma",              dForcingTermGamma);

     int iSolverOutput = 0;
     if (outputIters()) {
          if (uFlags & VERBOSE_MODE)           iSolverOutput |= NOX::Utils::Warning;
          if (uFlags & PRINT_CONVERGENCE_INFO) iSolverOutput |=
               NOX::Utils::OuterIteration | NOX::Utils::InnerIteration |
               NOX::Utils::LinearSolverDetails | NOX::Utils::StepperIteration |
               NOX::Utils::StepperDetails | NOX::Utils::Parameters |
               NOX::Utils::Details | NOX::Utils::OuterIterationStatusTest |
               NOX::Utils::TestDetails;
     }
     oPrintParam.set("Output Information", iSolverOutput);

     NOX::Abstract::PrePostOperator& oPrePost = *this;
     oSolverParam.sublist("Solver Options")
          .set("User Defined Pre/Post Operator", Teuchos::rcpFromRef(oPrePost));
     if (bInDerivativeSolver) {
          oSolverParam.sublist("Solver Options")
               .set("Status Test Check Type", "Complete");
     }

     static constexpr char szNLSolver[] = "Nonlinear Solver";

     if (uFlags & SOLVER_LINESEARCH_BASED) {
          oSolverParam.set(szNLSolver, "Line Search Based");
          auto& oLSP = oSolverParam.sublist("Line Search");
          std::string strLSMethod;
          if      (uFlags & LINESEARCH_BACKTRACK)    strLSMethod = "Backtrack";
          else if (uFlags & LINESEARCH_POLYNOMIAL)   strLSMethod = "Polynomial";
          else if (uFlags & LINESEARCH_MORE_THUENTE) strLSMethod = "More'-Thuente";
          else                                       strLSMethod = "Full Step";
          oLSP.set("Method", strLSMethod);
          if      (uFlags & SUFFICIENT_DEC_COND_ARMIJO_GOLDSTEIN)
               oLSP.sublist(strLSMethod).set("Sufficient Decrease Condition","Armijo-Goldstein");
          else if (uFlags & SUFFICIENT_DEC_COND_ARED_PRED)
               oLSP.sublist(strLSMethod).set("Sufficient Decrease Condition","Ared/Pred");
          auto& oLSM = oLSP.sublist(strLSMethod);
          oLSM.set("Max Iters",     iMaxIterLineSearch);
          oLSM.set("Minimum Step",  dMinStep);
          oLSM.set("Recovery Step", dRecoveryStep);
          if      (uFlags & RECOVERY_STEP_TYPE_CONST)
               oLSM.set("Recovery Step Type","Constant");
          else if (uFlags & RECOVERY_STEP_TYPE_LAST_STEP)
               oLSM.set("Recovery Step Type","Last Computed Step");
     } else if (uFlags & SOLVER_TRUST_REGION_BASED) {
          oSolverParam.set(szNLSolver, "Trust Region Based");
     } else if (uFlags & SOLVER_INEXACT_TRUST_REGION_BASED) {
          oSolverParam.set(szNLSolver, "Inexact Trust Region Based");
     } else if (uFlags & SOLVER_TENSOR_BASED) {
          oSolverParam.set(szNLSolver, "Tensor Based");
          oSolverParam.sublist("Line Search").set("Method","Curvilinear")
               .sublist("Curvilinear").set("Minimum Step", dMinStep);
     }

     if (uFlags & DIRECTION_NEWTON) {
          oDirectionParam.set("Method", "Newton");
          if      (uFlags & FORCING_TERM_CONSTANT) oNewtonParam.set("Forcing Term Method","Constant");
          else if (uFlags & FORCING_TERM_TYPE1)    oNewtonParam.set("Forcing Term Method","Type 1");
          else if (uFlags & FORCING_TERM_TYPE2)    oNewtonParam.set("Forcing Term Method","Type 2");
     } else if (uFlags & DIRECTION_STEEPEST_DESCENT) {
          oDirectionParam.set("Method", "Steepest Descent");
     } else if (uFlags & DIRECTION_NONLINEAR_CG) {
          oDirectionParam.set("Method", "NonlinearCG");
     } else if (uFlags & DIRECTION_BROYDEN) {
          oDirectionParam.set("Method", "Broyden");
          auto& oBroyden = oDirectionParam.sublist("Broyden");
          oBroyden.set("Restart Frequency", iIterationsBeforeAssembly);
          if      (uFlags & FORCING_TERM_CONSTANT) oBroyden.set("Forcing Term Method","Constant");
          else if (uFlags & FORCING_TERM_TYPE1)    oBroyden.set("Forcing Term Method","Type 1");
          else if (uFlags & FORCING_TERM_TYPE2)    oBroyden.set("Forcing Term Method","Type 2");
     }

     /* ---- NOX::Thyra::Group -------------------------------------------- */
     Teuchos::RCP<const Thyra::ModelEvaluator<TpetraSC>> pModelEvalConst = pModelEval;
     Teuchos::RCP<Thyra::LinearOpBase<TpetraSC>>         pJacobian = pModelEval->create_W_op();
     Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<TpetraSC>>
          pLOWSFactoryConst = pLOWSFactory;

     // Provide the assembled-Jacobian-based preconditioner.  For JFNK
     // this is essential (Belos GMRES needs it); for explicit-matrix +
     // iterative solver (Trust Region + GMRES), NOX may call
     // applyRightPreconditioning which also needs a preconditioner.
     Teuchos::RCP<Thyra::PreconditionerBase<TpetraSC>> pPrecOp =
          pModelEval->create_W_prec();
     Teuchos::RCP<Thyra::PreconditionerFactoryBase<TpetraSC>> pNullPrecFactory;

     auto grpPtr = Teuchos::rcp(new MBDynThyraGroup(
          *pSolutionView, pModelEvalConst, pJacobian,
          pLOWSFactoryConst, pPrecOp, pNullPrecFactory));

     /* ---- Status tests -------------------------------------------------- */
     auto converged =
          Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
     converged->addStatusTest(Teuchos::rcpFromRef(oResTest));
     if (oSolTest.dGetTolerance() > 0.) {
          converged->addStatusTest(Teuchos::rcpFromRef(oSolTest));
     }
     if (dWrmsRelTol > 0. && dWrmsAbsTol > 0. && !bInDerivativeSolver) {
          converged->addStatusTest(
               Teuchos::rcp(new NOX::StatusTest::NormWRMS(dWrmsRelTol, dWrmsAbsTol)));
     }

     auto pCombCriteria =
          Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
     pCombCriteria->addStatusTest(Teuchos::rcp(new NOX::StatusTest::FiniteValue));
     pCombCriteria->addStatusTest(converged);
     pCombCriteria->addStatusTest(
          Teuchos::rcp(new NOX::StatusTest::MaxIters(iMaxIter_a)));

     ForcePrecondRebuild();
     pNonlinearSolver = NOX::Solver::buildSolver(
          grpPtr, pCombCriteria, Teuchos::rcpFromRef(oSolverParam));
}

void NoxNonlinearSolver::OutputIteration(integer iIterCnt, bool bJacobian) const
{
#ifdef USE_MPI
     if (!bParallel || MBDynComm.Get_rank() == 0)
#endif
     {
          if (!(outputIters() || outputSolverConditionNumber())) return;
          silent_cout("\tIteration(" << iIterCnt << ") " << oResTest.dGetTest());
          if (bJacobian) {
               silent_cout(" J");
               if (outputSolverConditionNumber()) {
                    silent_cout(" cond=");
                    doublereal dCond;
                    if (pSolutionManager->bGetConditionNumber(dCond)) {
                         silent_cout(dCond);
                         if (outputSolverConditionStat()) {
                              AddCond(dCond);
                              silent_cout(" " << dGetCondMin()
                                         << " " << dGetCondMax()
                                         << " " << dGetCondAvg());
                         }
                    } else { silent_cout("NA"); }
               }
          }
          if (outputCPUTime()) {
               typedef std::chrono::duration<float,std::ratio<1,1>> FloatSec;
               auto flags = std::cout.flags(); auto prec = std::cout.precision();
               std::cout.setf(std::ios::scientific); std::cout.precision(2);
               silent_cout(" CPU:"
                    << FloatSec(dGetTimeCPU(CPU_RESIDUAL)).count()   << '+'
                    << FloatSec(dGetTimeCPU(CPU_JACOBIAN)).count()   << '+'
                    << FloatSec(dGetTimeCPU(CPU_LINEAR_SOLVER)).count());
               std::cout.flags(flags); std::cout.precision(prec);
          }
          if (oSolTest.getStatus() != NOX::StatusTest::Unevaluated) {
               silent_cout("\n\t\tSolErr " << oSolTest.dGetTest());
          }
          silent_cout('\n');
     }
}

void NoxNonlinearSolver::Attach(Solver* pS, const NonlinearProblem* pNLP)
{
     pSolver = pS;
     pSolutionManager = pSolver->pGetSolutionManager();
     if (pNLP != pNonlinearProblem) {
          DEBUGCERR("Resetting nonlinear solver\n");
          pNonlinearSolver.reset();
          bUpdateJacobian = true;
          ResetCond();
          if (pNonlinearProblem) bInDerivativeSolver = false;
     }
     if (!bKeepJacAcrossSteps) bUpdateJacobian = true;
     pNonlinearProblem = pNLP;
     VectorHandler* const pSol = pSolutionManager->pSolHdl();
     Size = pSol->iGetSize();
     pSol->Reset();
     DeltaX.ResizeReset(Size);
     XPrev.ResizeReset(Size);
     TmpRes.ResizeReset(Size);
}

void NoxNonlinearSolver::Residual(const VectorHandler* const pSol,
                                   VectorHandler* const pRes)
{
     CPUTimeGuard oCPUTimeRes(*this, CPU_RESIDUAL);
     // std::cerr << "DENTRO RESIDUAL" << std::endl;
     DeltaX.ScalarAddMul(*pSol, XPrev, -1.);
     XPrev = *pSol;
     pNonlinearProblem->Update(&DeltaX);
     pRes->Reset();
     VectorHandler* const pAbsRes = pGetResTest()->GetAbsRes();
     if (pAbsRes) pAbsRes->Reset();
     try {
          pNonlinearProblem->Residual(pRes, pAbsRes);
     } catch (const SolutionDataManager::ChangedEquationStructure&) {
          DEBUGCERR("Caught exception ChangedEquationStructure ...\n");
          if (bHonorJacRequest) {
               DEBUGCERR("Force update of preconditioner ...\n");
               ForcePrecondRebuild();
          }
     }
}

void NoxNonlinearSolver::Jacobian()
{
     CPUTimeGuard oCPUTimeJac(*this, CPU_JACOBIAN);
     bool bDone = false;
     MatrixHandler* pJac = nullptr;
     do {
          try {
               pJac = pSolutionManager->pMatHdl();
               pNonlinearProblem->Jacobian(pJac);
               pJac->PacMat();
               bDone = true;
          } catch (const MatrixHandler::ErrRebuildMatrix&) {
               silent_cout("NoxNonlinearSolver: rebuilding matrix...\n");
               pSolutionManager->MatrInitialize();
               // MBDynJacobianOp reads pMatHdl() at every apply() call, so
               // no pointer invalidation of the cached Jacobian op is needed
               // (unlike the old TpetraLinearOp path that cached a raw pointer).
          }
     } while (!bDone);
     ASSERT(pJac != nullptr);
#ifdef USE_MPI
     if (!bParallel || MBDynComm.Get_rank() == 0)
#endif
     {
          if (outputJac()) {
               silent_cout("Jacobian:\n");
               if (silent_out) pJac->Print(std::cout, MatrixHandler::MAT_PRINT_TRIPLET);
          }
     }
     ++TotJac;
}

void NoxNonlinearSolver::runPreLineSearch(const NOX::Solver::Generic&)
{
     DEBUGCERR("runPreLineSearch()\n");
     bInLineSearch = true;
     SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.);
     SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, 0);
}

void NoxNonlinearSolver::runPostLineSearch(const NOX::Solver::Generic&)
{
     DEBUGCERR("runPostLineSearch()\n");
     bInLineSearch = false;
     SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.);
     SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, 0);
}

} // anonymous namespace

/* =========================================================================
 * Public factory
 * ========================================================================= */
NonlinearSolver*
pAllocateNoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                            const NoxSolverParameters& oParam)
{
     NoxNonlinearSolver* pNLS = nullptr;
     SAFENEWWITHCONSTRUCTOR(pNLS, NoxNonlinearSolver,
                            NoxNonlinearSolver(oSolverOpt, oParam));
     return pNLS;
}

#endif /* USE_TRILINOS */
