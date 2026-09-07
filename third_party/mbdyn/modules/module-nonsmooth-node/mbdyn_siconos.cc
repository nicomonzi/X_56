/* $Header$ */
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
 * Author: Matteo Fancello <matteo.fancello@gmail.com>
 * Nonsmooth dynamics element;
 * uses SICONOS <http://siconos.gforge.inria.fr/>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <myassert.h>

#include "mbdyn_siconos.h"

#include <numerics/SolverOptions.h>
#include <numerics/LCP_Solvers.h>
#include <numerics/NumericsMatrix.h>
#include <numerics/NonSmoothDrivers.h>
#include <numerics/LinearComplementarityProblem.h>
#include <numerics/SiconosNumerics.h>

void
mbdyn_siconos_LCP_call(int sizep, double W_NN[], double bLCP[], double Pkp1[], double wlem[], solver_parameters& solparam)
{
        solparam.info = 0;

        LCPsolver lcpsol = solparam.solver;
        double tolerance = solparam.solvertol;
        int maxiternum = solparam.solveritermax;

        // LCP problem description:
        LinearComplementarityProblem OSNSProblem{};
        OSNSProblem.size = sizep;
        OSNSProblem.q = bLCP;

        NumericsMatrix MM{};

        MM.storageType = NM_DENSE;
        MM.matrix0 = W_NN;              // W_delassus
        MM.size0 = sizep;
        MM.size1 = sizep;
        OSNSProblem.M = &MM;

        SolverOptions* numerics_solver_options = solver_options_create(lcpsol);

        if (!numerics_solver_options) {
                silent_cerr("Unable to get default solver options for Siconos solver\n");
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        numerics_solver_options->iparam[0] = maxiternum;
        numerics_solver_options->dparam[0] = tolerance;

        solparam.info = linearComplementarity_driver(&OSNSProblem, Pkp1, wlem, numerics_solver_options);

        if (solparam.info) {
                silent_cerr("Siconos solver failed with status " << solparam.info << "\n");
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        solparam.processed_iterations = numerics_solver_options->iparam[1];
        solparam.resulting_error = numerics_solver_options->dparam[1];

        solver_options_delete(numerics_solver_options);
        free(numerics_solver_options);
}
