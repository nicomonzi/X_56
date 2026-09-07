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
	Copyright (C) 2024(-2025) all rights reserved.

	The copyright of this code is transferred
	to Pierangelo Masarati and Paolo Mantegazza
	for use in the software MBDyn as described
	in the GNU Public License version 2.1
*/

#include <cassert>
#include <iostream>
#include <iomanip>
#include <array>

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "myassert.h"
#include "mynewmem.h"
#include "solidcsl_impl.h"
#include "solidshape.h"
#include "demangle.h"

#ifdef USE_MPI
#include "mbcomm.h"
MPI::Intracomm MBDynComm = MPI::COMM_SELF; // FIXME: Workaround an undefined reference
#endif

template<typename function_t>
doublereal quad(const function_t& f, const doublereal a, const doublereal b, integer N)
{
     if ( a == b )
	  return 0.;

     N += N & 1;

     doublereal Y2n0 = 0., Y2n1 = 0.;

     for (integer i = 1; i <= N - 1; i += 2)
	  Y2n1 += f(a + (b - a) * i / N);

     for (integer i = 2; i <= N - 2; i += 2)
	  Y2n0 += f(a + (b - a) * i / N);

     return (b - a) / (3 * N) * (f(a) + f(b) + 4 * Y2n1 + 2 * Y2n0);
}

Mat3x3 CauchyStressTensor(const Mat3x3& F, ConstitutiveLaw6D& oCSL)
{
     const Mat3x3 C = F.Transpose() * F;
     const Mat3x3 G = (C - Eye3) / 2.;
     const Vec6 Eps(G(1, 1), G(2, 2), G(3, 3), 2. * G(1, 2), 2. * G(2, 3), 2. * G(3, 1));

     oCSL.Update(Eps);

     const Vec6& sigma = oCSL.GetF();

     Mat3x3 S;

     S(1, 1) = sigma(1);
     S(2, 2) = sigma(2);
     S(3, 3) = sigma(3);
     S(1, 2) = S(2, 1) = sigma(4);
     S(2, 3) = S(3, 2) = sigma(5);
     S(3, 1) = S(1, 3) = sigma(6);

     return F * S * F.Transpose() / F.dDet();
}

Mat3x3 CauchyStressTensor(const Mat3x3& F, ConstitutiveLaw9D& oCSL)
{
     const Mat3x3 G = (F.Transpose() * F - Eye3) / 2.;
     const Vec9 Eps{F(1, 1), F(2, 2), F(3, 3), F(1, 2), F(2, 1), F(1, 3), F(3, 1), F(2, 3), F(3, 2)};

     oCSL.Update(Eps);

     const Vec9& sigma = oCSL.GetF();

     Mat3x3 P;

     P(1, 1) = sigma(1);
     P(2, 2) = sigma(2);
     P(3, 3) = sigma(3);
     P(1, 2) = sigma(4);
     P(2, 1) = sigma(5);
     P(1, 3) = sigma(6);
     P(3, 1) = sigma(7);
     P(2, 3) = sigma(8);
     P(3, 2) = sigma(9);

     const Mat3x3 S = F.Inv() * P;

     return P * F.Transpose() / F.dDet();
}

template <typename ConstLawType>
Mat3x3 SedlanStress(ConstLawType& oCSL, doublereal r, doublereal D, doublereal lambda)
{
     Mat3x3 F = Zero3x3;

     F(1, 1) = F(2, 2) = 1. / sqrt(lambda);
     F(3, 3) = lambda;
     F(2, 3) = r * D;

     return CauchyStressTensor(F, oCSL);
}

template <typename ConstLawType>
doublereal SedlanTorque(ConstLawType& oCSL, doublereal r, doublereal D, doublereal lambda)
{
     const Mat3x3 T = SedlanStress(oCSL, r, D, lambda);

     const double Tphiz = T(2, 3);

     const double dm = 2. * M_PI * r * r * Tphiz;

     return dm;
}

template <typename ConstLawType>
doublereal SedlanForce(ConstLawType& oCSL, doublereal r, doublereal D, doublereal lambda)
{
     const Mat3x3 T = SedlanStress(oCSL, r, D, lambda);

     const double Tzz = T(3, 3);

     const double dn = 2. * M_PI * r * Tzz;

     return dn;
}

doublereal SedlanForceReference(doublereal C1, doublereal C2,  doublereal Ra, doublereal D, doublereal lambda)
{
     // Maxima expression:
     // N=(%pi*Ra^2*(4*C1*lambda^4+4*C2*lambda^3-C1*D^2*Ra^2*lambda-4*C1*lambda-2*C2*D^2*Ra^2-4*C2))/(3*lambda^3)$
     return M_PI * std::pow(Ra, 2) * (4. * C1 * std::pow(lambda, 4) + 4. * C2 * std::pow(lambda, 3) - C1 * std::pow(D, 2) * std::pow(Ra, 2) * lambda - 4. * C1 * lambda - 2. * C2 * std::pow(D, 2) * std::pow(Ra, 2) - 4. * C2) / (3. * std::pow(lambda, 3));
}

doublereal SedlanTorqueReference(doublereal C1, doublereal C2,  doublereal Ra, doublereal D, doublereal lambda)
{
     // Maxima expression:
     // M=(%pi*D*Ra^4*(C1*lambda+C2))/lambda^2$
     return M_PI * std::pow(Ra, 4) * (C1 * lambda + C2) * D / std::pow(lambda, 2);
}

doublereal NeoHookeanForceReference(doublereal mu, doublereal Ra, doublereal D, doublereal lambda)
{
     // Maxima expression:
     // N=(%pi*Ra^2*mu*(lambda-1)*(lambda+1))/lambda
     return (M_PI*std::pow(Ra,2)*mu*(lambda-1)*(lambda+1))/lambda;
}

doublereal NeoHookeanTorqueReference(doublereal mu, doublereal Ra, doublereal D, doublereal lambda)
{
     // Maxima expression:
     // M=(%pi*D*Ra^4*mu)/(2*lambda)$
     return (M_PI*D*std::pow(Ra,4)*mu)/(2*lambda);
}

template <typename ConstLawType>
struct IsotropicElasticityHelper;

template<>
struct IsotropicElasticityHelper<ConstitutiveLaw6D> {
     static Mat6x6 TangentOperator(doublereal E, doublereal nu) {
	  const doublereal a = nu / (1 - nu);
	  const doublereal b = (1 - 2 * nu) / (2 * (1 - nu));
	  const doublereal c = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
	  const doublereal d = a * c;
	  const doublereal e = b * c;

	  return Mat6x6(c, d, d, 0, 0, 0,
			d, c, d, 0, 0, 0,
			d, d, c, 0, 0, 0,
			0, 0, 0, e, 0, 0,
			0, 0, 0, 0, e, 0,
			0, 0, 0, 0, 0, e);
     }
};

template<>
struct IsotropicElasticityHelper<ConstitutiveLaw9D> {
     static Mat9x9 TangentOperator(doublereal E, doublereal nu) {
	  const doublereal c = (E*(nu-1))/((nu+1)*(2*nu-1));
	  const doublereal d = -(E*nu)/((nu+1)*(2*nu-1));
	  const doublereal e = E/(2*(nu+1));
	  const doublereal n = 0.;

	  return Mat9x9{c, d, d, n, n, n, n, n, n,
		    d, c, d, n, n, n, n, n, n,
		    d, d, c, n, n, n, n, n, n,
		    n, n, n, e, e, n, n, n, n,
		    n, n, n, e, e, n, n, n, n,
		    n, n, n, n, n, e, e, n, n,
		    n, n, n, n, n, e, e, n, n,
		    n, n, n, n, n, n, n, e, e,
		    n, n, n, n, n, n, n, e, e};
     }
};


template <typename ConstitutiveLawType>
void CheckConstitutiveLaw(ConstitutiveLawType& oMaterial, const typename ConstitutiveLawType::StressDerStrainType& Kref, const doublereal beta, const doublereal dTol)
{
     MBDYN_TESTSUITE_SCOPED_TRACE(oMaterial.GetName());

     typedef typename ConstitutiveLawType::StrainType StrainType;

     StrainType Eps = mb_zero<StrainType>(), EpsP = mb_zero<StrainType>();

     if constexpr(StrainType::iNumRowsStatic == 9) {
	       for (integer i = 1; i <= 3; ++i) {
		    Eps(i) = 1.;
	       }
     }

     const doublereal deltaEps = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.8);
     const doublereal deltaEpsP = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.8);

     std::cout << "\n" << oMaterial.GetName() << "\n";
     CheckTangentOperator(oMaterial, Eps, EpsP, Kref, beta, dTol, deltaEps, deltaEpsP);
}

template <typename TStress, typename TStressDerStrain, typename TStrain>
void SedlanConstLawTest(ConstitutiveLaw<TStress, TStressDerStrain, TStrain>& oCSL, const doublereal C1, const doublereal C2, const doublereal kappa)
{
     // Test case for hyperelastic material
     //
     // Konstantin Sedlan
     // Viskoelastisches Materialverhalten von Elastomerwerkstoffen
     // Experimentelle Untersuchung und Modellbildung
     // Kassel 2000
     //
     // Page 17, figure (3.1), page 18, equation (3.11) page 102, equation (4.114)
     //
     // https://kobra.uni-kassel.de/themes/Mirage2/scripts/mozilla-pdf.js/web/viewer.html?file=/bitstream/handle/123456789/2007110219579/DissertationKonstantinSedlan.pdf

     constexpr integer iNumSteps = 10;
     constexpr integer iNumInteg = 100;
     constexpr double Dmax = 1.;
     constexpr double lambdamin = 1.;
     constexpr double lambdamax = 10.;
     constexpr doublereal Ra = 25e-3;

     const doublereal delta = C2 / C1;
     const doublereal mu = C1 * (2. * (1. + delta));
     const doublereal E = (9 * kappa * mu) / (mu + 3 * kappa);
     const doublereal nu = (3 * kappa - 2 * mu) / (2 * mu + 6 * kappa);

     std::array<doublereal, iNumSteps + 1> N, M, Nref, Mref, lambda, D;

     for (integer i = 0; i <= iNumSteps; ++i) {
	  D[i] = i * Dmax / iNumSteps;
	  lambda[i] = i * (lambdamax - lambdamin) / iNumSteps + lambdamin;

	  const double ra = Ra / sqrt(lambda[i]);

	  auto dM = [&](const double r) {
	    return SedlanTorque(oCSL, r, D[i], lambda[i]);
	  };

	  auto dN = [&](const double r) {
			  return SedlanForce(oCSL, r, D[i], lambda[i]);
	  };

	  M[i] = quad(dM, 0., ra, iNumInteg);
	  N[i] = quad(dN, 0., ra, iNumInteg);

	  if (oCSL.GetName() == "NeoHookean6D") {
	       Nref[i] = NeoHookeanForceReference(mu, Ra, D[i], lambda[i]);
	       Mref[i] = NeoHookeanTorqueReference(mu, Ra, D[i], lambda[i]);
	  } else {
	       Nref[i] = SedlanForceReference(C1, C2, Ra, D[i], lambda[i]);
	       Mref[i] = SedlanTorqueReference(C1, C2, Ra, D[i], lambda[i]);
	  }
     }

     constexpr doublereal dTol = 1e-8;

     const auto flags = std::cout.flags();
     const auto precision = std::cout.precision();

     std::cout << oCSL.GetName() << ": actual values\n";

     std::cout << "\n-------------------------------------------------------------\n";
     std::cout << std::setw(14) << std::left << "D" << " |" << std::setw(14) << std::left << "lambda" << " |"
	       << std::setw(14) << std::left << "M" << " |" << std::setw(14) << std::left << "N" << "\n";
     std::cout << "-------------------------------------------------------------\n";

     for (integer i = 0; i <= iNumSteps; ++i) {
	  std::cout << std::setw(14) << std::right << std::fixed << std::setprecision(2) << D[i] << " |"
		    << std::setw(14) << std::right << std::fixed << std::setprecision(2) << lambda[i] << " |"
		    << std::fixed << std::right << std::setprecision(4) << std::setw(14) << M[i] << " |"
		    << std::fixed << std::right << std::setprecision(4) << std::setw(14) << N[i] << "\n";
     }

     std::cout << oCSL.GetName() << ": reference values\n";

     std::cout << "\n-------------------------------------------------------------\n";
     std::cout << std::setw(14) << std::left << "D" << " |" << std::setw(14) << std::left << "lambda" << " |"
	       << std::setw(14) << std::left << "M" << " |" << std::setw(14) << std::left << "N" << "\n";
     std::cout << "-------------------------------------------------------------\n";

     for (integer i = 0; i <= iNumSteps; ++i) {
	  std::cout << std::setw(14) << std::right << std::fixed << std::setprecision(2) << D[i] << " |"
		    << std::setw(14) << std::right << std::fixed << std::setprecision(2) << lambda[i] << " |"
		    << std::fixed << std::right << std::setprecision(4) << std::setw(14) << Mref[i] << " |"
		    << std::fixed << std::right << std::setprecision(4) << std::setw(14) << Nref[i] << "\n";
     }

     std::cout << "-------------------------------------------------------------\n";

     std::cout.flags(flags);
     std::cout.precision(precision);

     for (integer i = 1; i <= iNumSteps; ++i) {
	  MBDYN_TESTSUITE_ASSERT(std::fabs(M[i] / Mref[i] - 1.) < dTol);
	  MBDYN_TESTSUITE_ASSERT(std::fabs(N[i] / Nref[i] - 1.) < dTol);
     }

     typedef ConstitutiveLaw<TStress, TStressDerStrain, TStrain> ConstLawType_template;

     const auto Href = IsotropicElasticityHelper<ConstLawType_template>::TangentOperator(E, nu);

     const doublereal dTolTangentOperator = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.7) * E;

     CheckConstitutiveLaw(oCSL, Href, 0., dTolTangentOperator);
}

MBDYN_TESTSUITE_TEST(solidcsltest, SedlanTestMooneyRivlin6D)
{
     typedef ConstLawPreStress<Vec6, ConstLawPreStressType::CSL_PRE_NONE> PreStressNone;
     constexpr PreStressNone oPreStressNone;

     constexpr doublereal E = 100e6;
     constexpr doublereal nu = 0.3;
     constexpr doublereal delta = 0.25;
     constexpr doublereal mu = E / (2. * (1. + nu));
     constexpr doublereal kappa = E / (3. * (1. - 2. * nu));
     constexpr doublereal C1 = mu / (2. * (1. + delta));
     constexpr doublereal C2 = delta * C1;

     MooneyRivlinElastic<PreStressNone, PreStressNone> oMooneyRivlin6D(C1, C2, kappa, oPreStressNone, oPreStressNone);

     oMooneyRivlin6D.PutName("MooneyRivlin6D");

     SedlanConstLawTest(oMooneyRivlin6D, C1, C2, kappa);
}

MBDYN_TESTSUITE_TEST(solidcsltest, SedlanTestMooneyRivlin9D)
{
     typedef ConstLawPreStress<Vec6, ConstLawPreStressType::CSL_PRE_NONE> PreStressNone;
     constexpr PreStressNone oPreStressNone;

     constexpr doublereal E = 100e6;
     constexpr doublereal nu = 0.3;
     constexpr doublereal delta = 0.25;
     constexpr doublereal mu = E / (2. * (1. + nu));
     constexpr doublereal kappa = E / (3. * (1. - 2. * nu));
     constexpr doublereal C1 = mu / (2. * (1. + delta));
     constexpr doublereal C2 = delta * C1;

     MooneyRivlinElasticDefGrad<PreStressNone, PreStressNone> oMooneyRivlin9D(C1, C2, kappa, oPreStressNone, oPreStressNone);

     oMooneyRivlin9D.PutName("MooneyRivlin9D");

     SedlanConstLawTest(oMooneyRivlin9D, C1, C2, kappa);
}

#ifdef USE_MFRONT
MBDYN_TESTSUITE_TEST(solidcsltest, SedlanTestSignorini9D)
{
     constexpr doublereal E = 100e6;
     constexpr doublereal nu = 0.3;
     constexpr doublereal delta = 0.25;
     constexpr doublereal mu = E / (2. * (1. + nu));
     constexpr doublereal kappa = E / (3. * (1. - 2. * nu));
     constexpr doublereal C1 = mu / (2. * (1. + delta));
     constexpr doublereal C2 = delta * C1;

     static constexpr char szMFrontLibPath[] = "libHyperelasticityBehaviours-generic.so";

     mgis::behaviour::FiniteStrainBehaviourOptions oFiniteStrainOpt;

     oFiniteStrainOpt.stress_measure = mgis::behaviour::FiniteStrainBehaviourOptions::PK1;
     oFiniteStrainOpt.tangent_operator = mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF;

     mgis::behaviour::Behaviour oBehaviour = mgis::behaviour::load(oFiniteStrainOpt,
								   szMFrontLibPath,
								   "Signorini",
								   mgis::behaviour::Hypothesis::TRIDIMENSIONAL);

     MFrontGenericInterfaceCSL<Vec9, Mat9x9> oSignorini(oBehaviour, new NullDriveCaller);

     oSignorini.PutName("Signorini9D");
     oSignorini.setMaterialProperty("C10", C1);
     oSignorini.setMaterialProperty("C01", C2);
     oSignorini.setMaterialProperty("C20", 0.);
     oSignorini.setMaterialProperty("K", kappa);

     SedlanConstLawTest(oSignorini, C1, C2, kappa);
}
#endif

MBDYN_TESTSUITE_TEST(solidcsltest, SedlanTestNeoHookean6D)
{
     typedef ConstLawPreStress<Vec6, ConstLawPreStressType::CSL_PRE_NONE> PreStressNone;
     constexpr PreStressNone oPreStressNone;

     constexpr doublereal E = 100e6;
     constexpr doublereal nu = 0.3;
     constexpr doublereal delta = 0;
     constexpr doublereal mu = E / (2. * (1. + nu));
     constexpr doublereal kappa = E / (3. * (1. - 2. * nu));
     constexpr doublereal Lambda = E * nu / ((1. + nu) * (1. - 2. * nu));
     constexpr doublereal C1 = mu / (2. * (1. + delta));
     constexpr doublereal C2 = delta * C1;

     NeoHookeanElastic<PreStressNone, PreStressNone> oNeoHookean6D(mu, Lambda, oPreStressNone, oPreStressNone);

     oNeoHookean6D.PutName("NeoHookean6D");
     SedlanConstLawTest(oNeoHookean6D, C1, C2, kappa);
}

template <typename ConstitutiveLawType>
void CheckTangentOperator(ConstitutiveLawType& oMaterial, typename ConstitutiveLawType::StrainType& Eps, typename ConstitutiveLawType::StrainType& EpsP, const typename ConstitutiveLawType::StressDerStrainType& Kref, const doublereal beta, const doublereal dTol, const doublereal deltaEps, const doublereal deltaEpsP)
{
     oMaterial.Update(Eps, EpsP);

     auto F0 = oMaterial.GetF();
     auto K = oMaterial.GetFDE();
     auto D = oMaterial.GetFDEPrime();

     for (integer i = 1; i <= Kref.iGetNumRows(); ++i) {
	  for (integer j = 1; j <= Kref.iGetNumCols(); ++j) {
	       std::cout << std::setw(10) << std::setprecision(4) << Kref(i, j) << " ";
	  }
	  std::cout << "\n";
     }

     std::cout << "\n";

     for (integer i = 1; i <= K.iGetNumRows(); ++i) {
	  for (integer j = 1; j <= K.iGetNumCols(); ++j) {
	       std::cout << std::setw(10) << std::setprecision(4) << K(i, j) << " ";
	  }
	  std::cout << "\n";
     }

     for (integer i = 1; i <= K.iGetNumRows(); ++i) {
	  for (integer j = 1; j <= K.iGetNumCols(); ++j) {
	       MBDYN_TESTSUITE_ASSERT(fabs(K(i, j) - Kref(i, j)) < dTol);
	       MBDYN_TESTSUITE_ASSERT(fabs(D(i, j) - beta * Kref(i, j)) < dTol);
	  }
     }

     for (integer j = 1; j <= Eps.iGetNumRows(); ++j) {
	  const doublereal dEps0 = Eps(j);
	  Eps(j) = dEps0 + deltaEps;

	  oMaterial.Update(Eps, EpsP);

	  const auto F = oMaterial.GetF();

	  for (integer i = 1; i <= F.iGetNumRows(); ++i) {
	       MBDYN_TESTSUITE_ASSERT(std::fabs(F(i) - F0(i) - deltaEps * Kref(i, j)) < dTol);
	  }

	  Eps(j) = dEps0;
     }

     for (integer j = 1; j <= EpsP.iGetNumRows(); ++j) {
	  EpsP(j) = deltaEpsP;

	  oMaterial.Update(Eps, EpsP);

	  const auto F = oMaterial.GetF();

	  for (integer i = 1; i <= F.iGetNumRows(); ++i) {
	       MBDYN_TESTSUITE_ASSERT(std::fabs(F(i) - F0(i) - deltaEpsP * beta * Kref(i, j)) < dTol);
	  }

	  EpsP(j) = 0.;
     }

     oMaterial.Update(Eps, EpsP);
}

MBDYN_TESTSUITE_TEST(solidcsltest, SmallStrainTest)
{
     constexpr doublereal E = 210000e6;
     constexpr doublereal nu = 0.3;

     constexpr doublereal beta = 0.01;
     constexpr doublereal delta = 0.25;
     constexpr doublereal mu = E / (2. * (1. + nu));
     constexpr doublereal Lambda = E * nu / ((1. + nu) * (1. - 2. * nu));
     constexpr doublereal kappa = E / (3. * (1. - 2. * nu));
     constexpr doublereal C1 = mu / (2. * (1. + delta));
     constexpr doublereal C2 = delta * C1;
     constexpr doublereal ET = 2000e6;
     constexpr doublereal sigmavy = 2000e6;

     const doublereal dTol = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.7) * E;

     const Mat6x6 Kref = IsotropicElasticityHelper<ConstitutiveLaw6D>::TangentOperator(E, nu);

#ifdef USE_MFRONT
     static constexpr char szMFrontLibPath[] = "libPlasticityBehaviours-generic.so";

     mgis::behaviour::Behaviour oBehaviour = mgis::behaviour::load(szMFrontLibPath, "IsotropicLinearHardeningPlasticity", mgis::behaviour::Hypothesis::TRIDIMENSIONAL);

     MFrontGenericInterfaceCSL<Vec6, Mat6x6> oMFrontMaterial6D(oBehaviour, new NullDriveCaller);

     oMFrontMaterial6D.setMaterialProperty("YoungModulus", E);
     oMFrontMaterial6D.setMaterialProperty("PoissonRatio", nu);
     oMFrontMaterial6D.setMaterialProperty("YieldStrength", sigmavy);

     oMFrontMaterial6D.PutName("MFrontSmallStrain");
#endif

     typedef ConstLawPreStress<Vec6, ConstLawPreStressType::CSL_PRE_NONE> PreStressNone;
     constexpr PreStressNone oPreStressNone;

     MooneyRivlinElastic<PreStressNone, PreStressNone> oMooneyRivlin6D(C1, C2, kappa, oPreStressNone, oPreStressNone);

     oMooneyRivlin6D.PutName("MooneyRivlin6D");

     NeoHookeanElastic<PreStressNone, PreStressNone> oNeoHookean6D(mu, Lambda, oPreStressNone, oPreStressNone);

     oNeoHookean6D.PutName("NeoHookean6D");

     NeoHookeanViscoelastic<PreStressNone, PreStressNone> oNeoHookeanVisco6D(mu, Lambda, beta, oPreStressNone, oPreStressNone);

     oNeoHookeanVisco6D.PutName("NeoHookeanVisco6D");

     HookeanLinearElasticIsotropic<PreStressNone, PreStressNone> oHookean6D(mu, Lambda, oPreStressNone, oPreStressNone);

     oHookean6D.PutName("Hookean6D");

     HookeanLinearViscoelasticIsotropic<PreStressNone, PreStressNone> oHookeanVisco6D(mu, Lambda, beta, oPreStressNone, oPreStressNone);

     oHookeanVisco6D.PutName("HookeanVisco6D");

     BilinearIsotropicHardening<ConstitutiveLaw6D, ConstLawType::ELASTIC, PreStressNone> oBilinearIsotropicHardening6D(E, nu, ET, sigmavy, oPreStressNone);

     oBilinearIsotropicHardening6D.PutName("BilinearIsotropicHardening6D");

#ifdef USE_MFRONT
     CheckConstitutiveLaw(oMFrontMaterial6D, Kref, 0., dTol);
#endif
     CheckConstitutiveLaw(oMooneyRivlin6D, Kref, 0., dTol);
     CheckConstitutiveLaw(oNeoHookean6D, Kref, 0., dTol);
     CheckConstitutiveLaw(oNeoHookeanVisco6D, Kref, beta, dTol);
     CheckConstitutiveLaw(oHookean6D, Kref, 0., dTol);
     CheckConstitutiveLaw(oHookeanVisco6D, Kref, beta, dTol);
     CheckConstitutiveLaw(oBilinearIsotropicHardening6D, Kref, 0., dTol);
}

template <typename ElementType, sp_grad::index_type iDim>
bool bCheckShapeFunction()
{
     std::cout << "element type: \"" << mbdyn_demangle(typeid(ElementType)) << "\"\n";

     bool bRes = true;
     using namespace sp_grad;

     SpColVectorA<doublereal, iDim> r;
     SpColVectorA<doublereal, ElementType::iNumNodes> h;

     const doublereal dTol = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.9);

     for (index_type i = 1; i <= ElementType::iNumNodes; ++i) {
	  std::cout << "node: " << i << "\n";

	  ElementType::NodalPosition(i, r);
	  ElementType::ShapeFunction(r, h);

	  std::cout << "r = {" << r << "}\n";
	  std::cout << "h = {" << h << "}\n";

	  for (index_type j = 1; j <= ElementType::iNumNodes; ++j) {
	       if (fabs(h(j) - (i == j)) > dTol) {
		    bRes = false;
	       }
	  }
     }

     return bRes;
}

template <typename ElementType, sp_grad::index_type iDim>
bool bCheckShapeFunctionUPC()
{
     std::cout << "element type: \"" << ElementType::ElementName() << "\"\n";

     bool bRes = true;
     using namespace sp_grad;

     SpColVectorA<doublereal, iDim> r;
     SpColVectorA<doublereal, ElementType::ElemTypeDisplacement::iNumNodes> h;
     SpColVectorA<doublereal, ElementType::ElemTypePressure::iNumNodes> g;

     for (index_type i = 1; i <= ElementType::ElemTypeDisplacement::iNumNodes; ++i) {
	  std::cout << "node: " << i << "\n";

	  ElementType::ElemTypeDisplacement::NodalPosition(i, r);
	  ElementType::ElemTypeDisplacement::ShapeFunction(r, h);
	  ElementType::ElemTypePressure::ShapeFunction(r, g);

	  std::cout << "r = {" << r << "}\n";
	  std::cout << "h = {" << h << "}\n";
	  std::cout << "g = {" << g << "}\n";

	  for (index_type j = 1; j <= ElementType::ElemTypeDisplacement::iNumNodes; ++j) {
	       if (h(j) != (i == j)) {
		    bRes = false;
	       }

	       if (i <= ElementType::ElemTypePressure::iNumNodes && j <= ElementType::ElemTypePressure::iNumNodes) {
		    if (g(j) != (i == j)) {
			 bRes = false;
		    }
	       }
	  }
     }

     return bRes;
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionQ4)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Quadrangle4, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionQ8)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Quadrangle8, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionQ9)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Quadrangle9, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionQ8r)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Quadrangle8r, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT6h)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Triangle6h, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT10)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Triangle10, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH8u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron8u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH8p)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron8p, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH20u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron20u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH27u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron27u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH20upc)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunctionUPC<Hexahedron20upc, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH20ur)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron20ur, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH20upcr)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunctionUPC<Hexahedron20upcr, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionP6u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Pentahedron6u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionP15u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Pentahedron15u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionP18u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Pentahedron18u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionP15upc)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunctionUPC<Pentahedron15upc, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT4u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Tetrahedron4u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT10u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Tetrahedron10u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT10upc)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunctionUPC<Tetrahedron10upc, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT20u)
{
      MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Tetrahedron20u, 3>()));
}

MBDYN_DEFINE_OPERATOR_NEW_DELETE

int main(int argc, char* argv[])
{
     MBDYN_TESTSUITE_INIT(&argc, argv);

     return MBDYN_RUN_ALL_TESTS();
}
