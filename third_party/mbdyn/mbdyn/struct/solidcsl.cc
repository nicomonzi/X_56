/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include "solidcsl_impl.h"

void InitSolidCSL()
{
     // linear elasticity
     SetCL6D("hookean" "linear" "elastic" "isotropic", new HookeanLinearElasticIsotropicRead);
     SetCL7D("hookean" "linear" "elastic" "isotropic", new HookeanLinearElasticIsotropicIncompressibleRead);
     SetCL6D("hookean" "linear" "viscoelastic" "isotropic", new HookeanLinearViscoelasticIsotropicRead);
     // Need to add "hookean" because "linear elastic isotropic" is already in use

     // hyperelasticity
     SetCL6D("neo" "hookean" "elastic", new NeoHookeanReadElastic);
     SetCL6D("neo" "hookean" "viscoelastic", new NeoHookeanReadViscoelastic);
     SetCL6D("mooney" "rivlin" "elastic", new MooneyRivlinReadElastic<MooneyRivlinElastic, ConstitutiveLaw6D>);
     SetCL9D("mooney" "rivlin" "elastic", new MooneyRivlinReadElastic<MooneyRivlinElasticDefGrad, ConstitutiveLaw9D>);
     SetCL7D("mooney" "rivlin" "elastic", new MooneyRivlinReadElastic<MooneyRivlinElasticIncompressible, ConstitutiveLaw7D>);

     // plasticity
     SetCL6D("bilinear" "isotropic" "hardening", new BilinearIsotropicHardeningRead<ConstitutiveLaw6D, BilinearIsotropicHardening6D, ConstitutiveLawRead<Vec6, Mat6x6>>);
     SetCL7D("bilinear" "isotropic" "hardening", new BilinearIsotropicHardeningRead<ConstitutiveLaw7D, BilinearIsotropicHardening7D, ConstitutiveLawRead<Vec7, Mat7x7>>);

     // visco elasticity
#if defined(USE_LAPACK) && defined(HAVE_DGETRF) && defined(HAVE_DGETRS) && defined(HAVE_DGETRI)
     SetCL3D("linear" "viscoelastic" "maxwell" "1", new LinearViscoelasticMaxwell1Read<Vec3, Mat3x3>);
     SetCL6D("linear" "viscoelastic" "maxwell" "1", new LinearViscoelasticMaxwell1Read<Vec6, Mat6x6>);
     SetCL3D("linear" "viscoelastic" "maxwell" "n", new LinearViscoelasticMaxwellNRead<Vec3, Mat3x3>);
     SetCL6D("linear" "viscoelastic" "maxwell" "n", new LinearViscoelasticMaxwellNRead<Vec6, Mat6x6>);
#endif

#ifdef USE_MFRONT
     SetCL6D("mfront" "small" "strain", new MFrontGenericInterfaceCSLRead<Vec6, Mat6x6>);
     SetCL9D("mfront" "finite" "strain", new MFrontGenericInterfaceCSLRead<Vec9, Mat9x9>);
#endif

     SetDC6D("green" "lagrange" "strain", new GreenLagrangeStrainTplDCRead);
}
