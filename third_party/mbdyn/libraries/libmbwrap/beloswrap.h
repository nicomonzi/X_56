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

  Tpetra port: AztecOO → Belos, Amesos → Amesos2.
  This header replaces aztecoowrap.h.
*/

#ifndef ___BELOS_SOLUTION_MANAGER_H__INCLUDED__
#define ___BELOS_SOLUTION_MANAGER_H__INCLUDED__

#ifdef USE_TRILINOS

#include "solman.h"
#include "mbcomm.h"

/*
 * Allocate a Belos-based iterative solution manager (replaces AztecOO).
 *   uSolverFlags – preconditioner selection via LinSol::SOLVER_FLAGS_*
 */
SolutionManager*
pAllocateBelosSolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer iNLD,
     integer iMaxIter,
     doublereal dTolRes,
     integer iVerbose,
     unsigned uSolverFlags);

/*
 * Allocate an Amesos2-based direct solution manager (replaces Amesos).
 *   uSolverFlags – direct solver selection via LinSol::SOLVER_FLAGS_*
 */
SolutionManager*
pAllocateAmesos2SolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer iNLD,
     integer iVerbose,
     unsigned uSolverFlags);

/*
 * Finalize Tpetra/Kokkos runtime.  Call after all Trilinos objects
 * have been destroyed so that Kokkos finalize hooks run while static
 * objects (Belos MultiVecPool) are still alive.
 */
void mbdyn_trilinos_finalize();

#endif  /* USE_TRILINOS */
#endif  /* ___BELOS_SOLUTION_MANAGER_H__INCLUDED__ */
