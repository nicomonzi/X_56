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
  Copyright (C) 2018(-2019) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#ifndef __PASTIX_SOLUTION_MANAGER_H__INCLUDED__
#define __PASTIX_SOLUTION_MANAGER_H__INCLUDED__

#ifdef USE_PASTIX

#include <iostream>
#include <vector>
#include <complex>

#include "myassert.h"
#include "mynewmem.h"
#include "ls.h"
#include "solman.h"
#include "spmapmh.h"
#include "sp_gradient_spmh.h"

#define MPI_COMM_WORLD 0

#define save_HAVE_SCHED_SETAFFINITY HAVE_SCHED_SETAFFINITY
#define save_HAVE_CLOCK_GETTIME HAVE_CLOCK_GETTIME
#define save_HAVE_UNISTD_H HAVE_UNISTD_H
#define save_HAVE_GETOPT_LONG HAVE_GETOPT_LONG
#define save_HAVE_GETOPT_H HAVE_GETOPT_H
#define save_HAVE_STRING_H HAVE_STRING_H

#undef HAVE_SCHED_SETAFFINITY
#undef HAVE_CLOCK_GETTIME
#undef HAVE_UNISTD_H
#undef HAVE_GETOPT_LONG
#undef HAVE_GETOPT_H
#undef HAVE_STRING_H

extern "C" {
#include <pastix.h>
}

#undef HAVE_SCHED_SETAFFINITY
#undef HAVE_CLOCK_GETTIME
#undef HAVE_UNISTD_H
#undef HAVE_GETOPT_LONG
#undef HAVE_GETOPT_H
#undef HAVE_STRING_H

#define HAVE_SCHED_SETAFFINITY save_HAVE_SCHED_SETAFFINITY
#define HAVE_CLOCK_GETTIME save_HAVE_CLOCK_GETTIME
#define HAVE_UNISTD_H save_HAVE_UNISTD_H
#define HAVE_GETOPT_LONG save_HAVE_GETOPT_LONG
#define HAVE_GETOPT_H save_HAVE_GETOPT_H
#define HAVE_STRING_H save_HAVE_STRING_H

class PastixSolver: public LinearSolver {
private:
     struct SpMatrix: spmatrix_t {
	  SpMatrix();
	  SpMatrix(const SpMatrix&)=delete;	 
	  ~SpMatrix();

	  SpMatrix& operator=(const SpMatrix&)=delete;

	  bool MakeCompactForm(const SparseMatrixHandler& mh);
	  
	  doublereal* pAx() const { return reinterpret_cast<doublereal*>(values); }
	  pastix_int_t* pAi() const { return rowptr; }
	  pastix_int_t* pAp() const { return colptr; }
	  integer Nz() const { return iNumNonZeros; }
     private:
	  template <typename T>
	  static inline T* pAllocate(T* pMem, size_t nSize);
	  void Allocate(size_t iNumNZ, size_t iMatSize);
	  integer iNumNonZeros;
     };
     
     pastix_int_t iparm[IPARM_SIZE];
     doublereal dparm[DPARM_SIZE];
     pastix_data_t* pastix_data;
     mutable SpMatrix spm;
     mutable bool bDoOrdering;
     
public:
     explicit PastixSolver(SolutionManager* pSM,
			   integer iDim,
			   integer iNumIter,
                           doublereal dTolRefine,
			   integer iNumThreads,
			   unsigned uSolverFlags,
			   doublereal dCompressTol,
			   doublereal dMinRatio,
			   integer iVerbose);
     ~PastixSolver();

#ifdef DEBUG
     void IsValid(void) const;
#endif /* DEBUG */

     void Initialize() { bDoOrdering = true; }
     virtual void Solve(void) const;
     SpMatrix& PastixMakeCompactForm(SparseMatrixHandler& mh);
};

template <typename MatrixHandlerType>
class PastixSolutionManager: public SolutionManager {
private:
    std::vector<doublereal> x;
    std::vector<doublereal> b;
    mutable MyVectorHandler xVH, bVH;
    ScaleOpt scale;
    MatrixScaleBase* pMatScale;

protected:
    mutable MatrixHandlerType A;

    PastixSolver* pGetSolver() { return static_cast<PastixSolver*>(pLS); }

    template <typename MH>
    void ScaleMatrixAndRightHandSide(MH &mh);

    template <typename MH>
    MatrixScale<MH>& GetMatrixScale();

    void ScaleSolution(void);

public:
    PastixSolutionManager(integer iDim,
			  integer iNumThreads,
			  integer iNumIter,
                          doublereal dTolRefine,
			  const ScaleOpt& scale = ScaleOpt(),
			  unsigned uSolverFlags = 0u,
			  doublereal dCompressTol = 0.,
			  doublereal dMinRatio = 1.,
			  integer iVerbose = 0);
    virtual ~PastixSolutionManager(void);
#ifdef DEBUG
    virtual void IsValid(void) const;
#endif /* DEBUG */
    virtual void MatrReset(void);
    virtual void MatrInitialize(void);
    virtual void Solve(void);
    virtual void MakeCompressedColumnForm(void);
    virtual MatrixHandler* pMatHdl(void) const;
    virtual VectorHandler* pResHdl(void) const;
    virtual VectorHandler* pSolHdl(void) const;
};

#endif

#endif
