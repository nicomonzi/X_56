/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
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
  */

#ifndef __NOX_SOLVER_H__INCLUDED__
#define __NOX_SOLVER_H__INCLUDED__

#ifdef USE_TRILINOS
#include <map>
#include <string>
#include <variant>
#include "nonlin.h"

struct NoxSolverParameters: public CommonNonlinearSolverParam {
     NoxSolverParameters();

     enum SolverFlags {
          JACOBIAN_NEWTON_KRYLOV = 0x10,
          JACOBIAN_NEWTON        = 0x20,
          JACOBIAN_OPERATOR_MASK = JACOBIAN_NEWTON_KRYLOV | JACOBIAN_NEWTON,
          SOLVER_LINESEARCH_BASED   = 0x40,
          SOLVER_TRUST_REGION_BASED = 0x80,
          SOLVER_INEXACT_TRUST_REGION_BASED = 0x100,
          SOLVER_TENSOR_BASED               = 0x200,
          SOLVER_MASK = SOLVER_LINESEARCH_BASED |
          SOLVER_TRUST_REGION_BASED |
          SOLVER_INEXACT_TRUST_REGION_BASED |
          SOLVER_TENSOR_BASED,
          DIRECTION_NEWTON = 0x400,
          DIRECTION_STEEPEST_DESCENT = 0x800,
          DIRECTION_NONLINEAR_CG = 0x1000,
          DIRECTION_BROYDEN = 0x2000,
          DIRECTION_MASK = DIRECTION_NEWTON |
          DIRECTION_STEEPEST_DESCENT |
          DIRECTION_BROYDEN,
          FORCING_TERM_CONSTANT =  0x4000,
          FORCING_TERM_TYPE1    =  0x8000,
          FORCING_TERM_TYPE2    = 0x10000,
          FORCING_TERM_MASK = FORCING_TERM_CONSTANT |
          FORCING_TERM_TYPE1 |
          FORCING_TERM_TYPE2,
          LINESEARCH_BACKTRACK    = 0x20000,
          LINESEARCH_POLYNOMIAL   = 0x40000,
          LINESEARCH_MORE_THUENTE = 0x80000,
          LINESEARCH_MASK = LINESEARCH_BACKTRACK |
          LINESEARCH_POLYNOMIAL |
          LINESEARCH_MORE_THUENTE,
          LINEAR_SOLVER_BLOCK_GMRES                       = 0x100000,
          LINEAR_SOLVER_PSEUDO_BLOCK_GMRES                = 0x200000,
          LINEAR_SOLVER_BLOCK_CG                          = 0x400000,
          LINEAR_SOLVER_PSEUDO_BLOCK_CG                   = 0x800000,
          LINEAR_SOLVER_BLOCK_STOCHASTIC_CG               = 0x1000000,
          LINEAR_SOLVER_GCRODR                            = 0x2000000,
          LINEAR_SOLVER_RCG                               = 0x4000000,
          LINEAR_SOLVER_MINRES                            = 0x8000000,
          LINEAR_SOLVER_TFQMR                             = 0x10000000,
          LINEAR_SOLVER_BICGSTAB                          = 0x20000000,
          LINEAR_SOLVER_FIXED_POINT                       = 0x40000000,
          LINEAR_SOLVER_TPETRA_GMRES                      = 0x80000000,
          LINEAR_SOLVER_TPETRA_GMRES_PIPELINE             = 0x100000000,
          LINEAR_SOLVER_TPETRA_GMRES_SINGLE_REDUCE        = 0x200000000,
          LINEAR_SOLVER_TPETRA_GMRES_SSTEP                = 0x400000000,
          LINEAR_SOLVER_MASK = LINEAR_SOLVER_BLOCK_GMRES |
          LINEAR_SOLVER_PSEUDO_BLOCK_GMRES |
          LINEAR_SOLVER_BLOCK_CG |
          LINEAR_SOLVER_PSEUDO_BLOCK_CG |
          LINEAR_SOLVER_BLOCK_STOCHASTIC_CG |
          LINEAR_SOLVER_GCRODR |
          LINEAR_SOLVER_RCG |
          LINEAR_SOLVER_MINRES |
          LINEAR_SOLVER_TFQMR |
          LINEAR_SOLVER_BICGSTAB |
          LINEAR_SOLVER_FIXED_POINT |
          LINEAR_SOLVER_TPETRA_GMRES |
          LINEAR_SOLVER_TPETRA_GMRES_PIPELINE |
          LINEAR_SOLVER_TPETRA_GMRES_SINGLE_REDUCE |
          LINEAR_SOLVER_TPETRA_GMRES_SSTEP,
          RECOVERY_STEP_TYPE_CONST     = 0x800000000,
          RECOVERY_STEP_TYPE_LAST_STEP = 0x1000000000,
          RECOVERY_STEP_TYPE_MASK = RECOVERY_STEP_TYPE_CONST |
          RECOVERY_STEP_TYPE_LAST_STEP,
          USE_PRECOND_AS_SOLVER = 0x2000000000,
          SUFFICIENT_DEC_COND_ARMIJO_GOLDSTEIN = 0x4000000000,
          SUFFICIENT_DEC_COND_ARED_PRED = 0x80000000000,
          SUFFICIENT_DEC_COND_MASK = SUFFICIENT_DEC_COND_ARMIJO_GOLDSTEIN |
                                      SUFFICIENT_DEC_COND_ARED_PRED
     };

     doublereal dWrmsRelTol;
     doublereal dWrmsAbsTol;
     doublereal dTolLinSol;
     doublereal dMinStep;
     doublereal dRecoveryStep;
     doublereal dForcingTermMinTol;
     doublereal dForcingTermMaxTol;
     doublereal dForcingTermAlpha;
     doublereal dForcingTermGamma;
     integer iMaxIterLinSol;
     integer iKrylovSubSpaceSize;
     integer iMaxIterLineSearch;
     integer iInnerIterBeforeAssembly;

     // -----------------------------------------------------------------------
     // Generic Belos solver-parameter overrides.
     //
     // At parse time, the "belos parameters" keyword opens a block where the
     // user can set any parameter that the chosen Belos solver accepts.  The
     // parameters are stored here as a map keyed on the Belos parameter name
     // (e.g. "Flexible Gmres", "Num Recycled Blocks", "Orthogonalization").
     // Values are typed variants matching the four scalar types that Belos
     // parameter lists use: bool, int, double, std::string.  Teuchos-internal
     // types (RCP<ostream> etc.) are intentionally excluded.
     //
     // The map is applied in noxsolver.cc::BuildSolver(), where it is
     // validated against the solver's getValidParameters() list before being
     // forwarded to Belos.  Parameters that collide with the named MBDyn
     // keywords (Maximum Iterations, Convergence Tolerance, Num Blocks,
     // Output Frequency, Verbosity, Output Style) are silently overridden by
     // the named keywords and should not be set here.
     // -----------------------------------------------------------------------
     using BelosParamValue = std::variant<bool, int, double, std::string>;
     std::map<std::string, BelosParamValue> oBelosParams;
};

NonlinearSolver*
pAllocateNoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                            const NoxSolverParameters& oParam);

// -------------------------------------------------------------------------
// BelosParamType – the scalar type of a Belos ParameterList entry.
//
// Used by the input parser (solver.cc) to dispatch to the correct
// HighParser getter without needing to include any Trilinos header.
// UNSETTABLE is returned for entries whose C++ type is not one of the four
// plain scalar types (e.g. Teuchos::RCP<std::ostream>, Teuchos::Array<…>).
// UNKNOWN is returned when the parameter name does not exist in the solver's
// valid-parameter list.
// -------------------------------------------------------------------------
enum class BelosParamType { BOOL, INT, DOUBLE, STRING, UNSETTABLE, UNKNOWN };

// Query the type of Belos parameter `sParamName` for solver `sSolverType`.
// This calls Belos::SolverFactory::create() with null params the first time
// each solver type is queried and caches the result; subsequent queries are
// O(log n) map lookups.
// Implementation: noxsolver.cc (requires Belos headers).
BelosParamType eGetBelosParamType(const std::string& sSolverType,
                                   const std::string& sParamName);

#endif
#endif
