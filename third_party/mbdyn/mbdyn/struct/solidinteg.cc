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

#include "solidinteg.h"

constexpr doublereal fabs_ce(doublereal x) {
        return x >= 0. ? x : -x;
}

namespace {
     constexpr doublereal dCollocationFunction1D(const doublereal x) {
          return x + 3. * x * x + 2.;
     }

     constexpr doublereal dCollocationFunction1DRes() {
          // maxima: integrate(1*x+3*x^2+2,x,-1,1);
          return 6.;
     }

     constexpr doublereal dCollocationTest1DRec(const doublereal r[], const doublereal alpha[], const int N) {
          return dCollocationFunction1D(r[0]) * alpha[0]
               + ((N > 1) ? dCollocationTest1DRec(r + 1, alpha + 1, N - 1) : 0.); // allow use to use -std=c++11
     }

     constexpr bool bCollocationTest1D(const doublereal r[], const doublereal alpha[], const int N, const doublereal dTol = 0.) {
          return fabs_ce(dCollocationTest1DRec(r, alpha, N) / dCollocationFunction1DRes() - 1.) <= dTol;
     }

     constexpr doublereal dCollocationFunction2D(const doublereal x, const doublereal y) {
          return 1. + x + y + x * y + x * x + y * y;
     }

     constexpr doublereal dCollocationFunction2DRes() {
          // maxima: integrate(subst([yi=y], integrate(1 + x + y + x * y + x * x + y * y, x, 0, 1 - yi)), y, 0, 1);
          return 25./24.;
     }

     constexpr doublereal dCollocationTest2DRec(const doublereal r[], const doublereal s[], const doublereal alpha[], const int N) {
          return dCollocationFunction2D(r[0], s[0]) * alpha[0]
               + ((N > 1) ? dCollocationTest2DRec(r + 1, s + 1, alpha + 1, N - 1) : 0.); // allow us to use -std=c++11
     }

     constexpr bool bCollocationTest2D(const doublereal r[], const doublereal s[], const doublereal alpha[], const int N, const doublereal dTol = 0.) {
          return fabs_ce(dCollocationTest2DRec(r, s, alpha, N) / dCollocationFunction2DRes() - 1.) <= dTol;
     }

     constexpr doublereal dCollocationFunction3D(const doublereal x, const doublereal y, const double z) {
          return 1. + x + y + z + x * x + y * y + z * z;
     }

     constexpr doublereal dCollocationTest3DRec(const doublereal r[], const doublereal s[], const doublereal t[], const doublereal alpha[], const int N) {
          return dCollocationFunction3D(r[0], s[0], t[0]) * alpha[0]
               + ((N > 1) ? dCollocationTest3DRec(r + 1, s + 1, t + 1, alpha + 1, N - 1) : 0.); // allow us to use -std=c++11
     }

     constexpr doublereal dCollocationFunction3DResTet() {
          // maxima: integrate(subst([zi=z], integrate(subst([yi=y], integrate(1. + x + y + z + x * x + y * y + z * z, x, 0, 1 - yi - zi)), y, 0, 1 - zi)), z, 0, 1);
          return 41. / 120.;
     }

     constexpr doublereal dCollocationFunction3DResPenta() {
          // maxima: integrate(subst([zi=z], integrate(subst([yi=y], integrate(1. + x + y + z + x * x + y * y + z * z, x, -1, 1)), y, 0, 1 - zi)), z, 0, 1);
          return 7. / 3.;
     }

     constexpr bool bCollocationTest3D(const doublereal r[], const doublereal s[], const doublereal t[], const doublereal alpha[], const int N, const double dRes, const doublereal dTol = 0.) {
          return fabs_ce(dCollocationTest3DRec(r, s, t, alpha, N) / dRes - 1.) <= dTol;
     }

     static_assert(bCollocationTest1D(Gauss2_1D::ri, Gauss2_1D::alphai, Gauss2_1D::iGaussOrder), "unit test for collocation rule failed");
     static_assert(bCollocationTest1D(Gauss3_1D::ri, Gauss3_1D::alphai, Gauss3_1D::iGaussOrder), "unit test for collocation rule failed");
     static_assert(bCollocationTest2D(CollocTria6h::zeta, CollocTria6h::eta, CollocTria6h::w, CollocTria6h::iNumEvalPoints), "unit test for collocation rule failed");
     static_assert(bCollocationTest2D(CollocTria10::r1, CollocTria10::r2, CollocTria10::w, CollocTria10::iNumEvalPoints, 5 * std::numeric_limits<doublereal>::epsilon()), "unit test for collocation rule failed");
     static_assert(bCollocationTest3D(CollocTet10h::r1, CollocTet10h::s1, CollocTet10h::t1, CollocTet10h::w1, CollocTet10h::N1, dCollocationFunction3DResTet()), "unit test for collocation rule failed");
     static_assert(bCollocationTest3D(CollocTet10h::r2, CollocTet10h::s2, CollocTet10h::t2, CollocTet10h::w2, CollocTet10h::N2, dCollocationFunction3DResTet(), std::numeric_limits<doublereal>::epsilon()), "unit test for collocation rule failed");
     static_assert(bCollocationTest3D(CollocTet10h::r3, CollocTet10h::s3, CollocTet10h::t3, CollocTet10h::w3, CollocTet10h::N3, dCollocationFunction3DResTet(), 5 * std::numeric_limits<doublereal>::epsilon()), "unit test for collocation rule failed");
     static_assert(bCollocationTest3D(CollocTet20::ri1, CollocTet20::si1, CollocTet20::ti1, CollocTet20::wi1, CollocTet20::N1, dCollocationFunction3DResTet(), std::numeric_limits<doublereal>::epsilon()), "unit test for collocation rule failed");
     static_assert(bCollocationTest3D(CollocTet20::ri3, CollocTet20::si3, CollocTet20::ti3, CollocTet20::wi3, CollocTet20::N3, dCollocationFunction3DResTet(), 15 * std::numeric_limits<doublereal>::epsilon()), "unit test for collocation rule failed");
     static_assert(bCollocationTest3D(CollocPenta18::ri, CollocPenta18::si, CollocPenta18::ti, CollocPenta18::wi, CollocPenta18::iNumEvalPointsStiffness, dCollocationFunction3DResPenta(), std::numeric_limits<doublereal>::epsilon()), "unit test for collocation rule failed");
}

constexpr sp_grad::index_type Gauss2_1D::iGaussOrder;
constexpr doublereal Gauss2_1D::ri[];
constexpr doublereal Gauss2_1D::alphai[];
constexpr doublereal Gauss2Lumped_1D::ri[];
constexpr doublereal Gauss2Lumped_1D::alphai[];

constexpr sp_grad::index_type IntegLayout2_2D::ridx[];
constexpr sp_grad::index_type IntegLayout2_2D::sidx[];

constexpr sp_grad::index_type IntegLayout3_2D::ridx[];
constexpr sp_grad::index_type IntegLayout3_2D::sidx[];

constexpr sp_grad::index_type IntegLayout2_3D::ridx[];
constexpr sp_grad::index_type IntegLayout2_3D::sidx[];
constexpr sp_grad::index_type IntegLayout2_3D::tidx[];

constexpr sp_grad::index_type IntegLayout3_3D::ridx[];
constexpr sp_grad::index_type IntegLayout3_3D::sidx[];
constexpr sp_grad::index_type IntegLayout3_3D::tidx[];

constexpr sp_grad::index_type CollocTria6h::iNumEvalPoints;
constexpr doublereal CollocTria6h::A;
constexpr doublereal CollocTria6h::B;
constexpr doublereal CollocTria6h::P1;
constexpr doublereal CollocTria6h::P2;
constexpr doublereal CollocTria6h::zeta[7];
constexpr doublereal CollocTria6h::eta[7];
constexpr doublereal CollocTria6h::w[7];

constexpr doublereal Gauss3_1D::ri[];
constexpr doublereal Gauss3_1D::alphai[];
constexpr doublereal Gauss3Lumped_1D::ri[];
constexpr doublereal Gauss3Lumped_1D::alphai[];

constexpr sp_grad::index_type CollocPenta15::M;
constexpr sp_grad::index_type CollocPenta15::N;
constexpr sp_grad::index_type CollocPenta15::M_lumped;
constexpr sp_grad::index_type CollocPenta15::iNumEvalPointsStiffness;
constexpr sp_grad::index_type CollocPenta15::iNumEvalPointsMass;
constexpr sp_grad::index_type CollocPenta15::iNumEvalPointsMassLumped;

constexpr doublereal CollocPenta15::ri[];
constexpr doublereal CollocPenta15::si[];
constexpr doublereal CollocPenta15::wi[];
constexpr doublereal CollocPenta15::ti[];
constexpr doublereal CollocPenta15::alphai[];

constexpr doublereal CollocPenta15::ri_lumped[];
constexpr doublereal CollocPenta15::si_lumped[];
constexpr doublereal CollocPenta15::wi_lumped[];
constexpr doublereal CollocPenta15::ti_lumped[];
constexpr doublereal CollocPenta15::alphai_lumped[];

constexpr sp_grad::index_type CollocPenta18::iNumEvalPointsStiffness;
constexpr sp_grad::index_type CollocPenta18::iNumEvalPointsMass;
constexpr sp_grad::index_type CollocPenta18::iNumEvalPointsMassLumped;

constexpr doublereal CollocPenta18::ri[];
constexpr doublereal CollocPenta18::si[];
constexpr doublereal CollocPenta18::ti[];
constexpr doublereal CollocPenta18::wi[];

constexpr doublereal CollocPenta18::ri_lumped[];
constexpr doublereal CollocPenta18::si_lumped[];
constexpr doublereal CollocPenta18::ti_lumped[];

constexpr sp_grad::index_type CollocTet10h::iNumEvalPointsStiffness;
constexpr sp_grad::index_type CollocTet10h::iNumEvalPointsMass;
constexpr sp_grad::index_type CollocTet10h::iNumEvalPointsMassLumped;

constexpr doublereal CollocTet10h::r1[];
constexpr doublereal CollocTet10h::s1[];
constexpr doublereal CollocTet10h::t1[];
constexpr doublereal CollocTet10h::w1[];

constexpr doublereal CollocTet10h::r2[];
constexpr doublereal CollocTet10h::s2[];
constexpr doublereal CollocTet10h::t2[];
constexpr doublereal CollocTet10h::w2[];

constexpr doublereal CollocTet10h::r3[];
constexpr doublereal CollocTet10h::s3[];
constexpr doublereal CollocTet10h::t3[];
constexpr doublereal CollocTet10h::w3[];
