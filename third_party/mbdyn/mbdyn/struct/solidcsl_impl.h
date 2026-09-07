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

#ifndef __SOLIDCSL_IMPL_H__INCLUDED__
#define __SOLIDCSL_IMPL_H__INCLUDED__
#include <memory>
#include "constltp.h"
#include "dataman.h"
#include "solidcsl.h"
#include "drive_.h"
#include "tensorformat.h"
#include <ac/lapack.h>

#ifdef USE_MFRONT
#include <constexpr_math.h>
// MFrontGenericInterfaceSupport by Thomas Helfer (https://thelfer.github.io/mgis/web/index.html)
// To be used with the MFront code generation tool (https://thelfer.github.io/tfel/web/index.html)
#include <MGIS/Behaviour/State.hxx>
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#endif

class GreenLagrangeStrainTplDC : public TplDriveCaller<Vec6> {
public:
     GreenLagrangeStrainTplDC(const DriveHandler* pDriveHdl, std::unique_ptr<TplDriveCaller<Vec6>>&& pStrainTmp)
	  :TplDriveCaller<Vec6>(pDriveHdl), pStrain(std::move(pStrainTmp)) {
	  DEBUGCOUTFNAME("GreenLagrangeStrainTplDC::GreenLagrangeStrainTplDC()\n");
     }

     virtual ~GreenLagrangeStrainTplDC() override {
	  DEBUGCOUTFNAME("GreenLagrangeStrainTplDC::~GreenLagrangeStrainTplDC()\n");
     }

     virtual TplDriveCaller<Vec6>* pCopy() const override {
	  using namespace sp_grad;

	  GreenLagrangeStrainTplDC* pDC = nullptr;

	  std::unique_ptr<TplDriveCaller<Vec6>> pStrainTmp(pStrain->pCopy());

	  SAFENEWWITHCONSTRUCTOR(pDC,
				 GreenLagrangeStrainTplDC,
				 GreenLagrangeStrainTplDC(pDrvHdl, std::move(pStrainTmp)));

	  return pDC;
     }

     virtual std::ostream& Restart(std::ostream& out) const override {
	  out << "green lagrange strain, ";
	  return Restart_int(out);
     }

     virtual std::ostream& Restart_int(std::ostream& out) const override {
	  return pStrain->Restart(out);
     }

     virtual Vec6 Get(const doublereal& dTime) const override {
	  using namespace sp_grad;
	  using TFS = TensorFormat::Symmetric3x3;

	  const Vec6 epsilonBar = pStrain->Get(dTime);

	  Vec6 epsilon;

	  for (index_type i = 1; i <= 3; ++i) {
	       epsilon(i) = 0.5 * (std::pow(epsilonBar(i) + 1., 2) - 1.);
	  }

	  for (index_type i = 1; i <= 3; ++i) {
	       epsilon(i + 3) = epsilonBar(i + 3) * (1. + epsilonBar(TFS::i1[i - 1 + 3])) * (1. + epsilonBar(TFS::i2[i - 1 + 3]));
	  }

	  return epsilon;
     }

     inline bool bIsDifferentiable() const override {
	  return pStrain->bIsDifferentiable();
     }

     virtual Vec6 GetP(const doublereal& dTime) const override {
	  using namespace sp_grad;
	  using TFS = TensorFormat::Symmetric3x3;

	  const Vec6 epsilonBar = pStrain->Get(dTime);
	  const Vec6 epsilonBarP = pStrain->GetP(dTime);

	  Vec6 epsilonP;

	  for (index_type i = 1; i <= 3; ++i) {
	       epsilonP(i) = (epsilonBar(i) + 1.) * epsilonBarP(i);
	  }

	  for (index_type i = 4; i <= 6; ++i) {
	       epsilonP(i) = epsilonBarP(i) * (1. + epsilonBar(TFS::i1[i - 1])) * (1. + epsilonBar(TFS::i2[i - 1]))
		    + epsilonBar(i) * (epsilonBarP(TFS::i1[i - 1]) * (1. + epsilonBar(TFS::i2[i - 1]))
				       + (1 + epsilonBar(TFS::i1[i - 1])) * epsilonBarP(TFS::i2[i - 1]));
	  }

	  return epsilonP;
     }

     inline int getNDrives() const override {
	  return pStrain->getNDrives();
     }

private:
     const std::unique_ptr<TplDriveCaller<Vec6>> pStrain;
};

struct GreenLagrangeStrainTplDCRead: TplDriveCallerRead<Vec6> {
	virtual TplDriveCaller<Vec6>*
	Read(const DataManager* pDM, MBDynParser& HP) override {
	     using namespace sp_grad;

	     std::unique_ptr<TplDriveCaller<Vec6>> pStrainTmp(HP.GetTplDriveCaller<Vec6>());

	     TplDriveCaller<Vec6>* pTplDC = nullptr;

	     SAFENEWWITHCONSTRUCTOR(pTplDC,
				    GreenLagrangeStrainTplDC,
				    GreenLagrangeStrainTplDC(pDM->pGetDrvHdl(), std::move(pStrainTmp)));

	     return pTplDC;
	}
};

enum class ConstLawPreStressType: integer {
     CSL_PRE_NONE,
     CSL_PRE_CONST,
     CSL_PRE_DRIVE
};

template <typename T, ConstLawPreStressType>
class ConstLawPreStress;

template <typename T>
class ConstLawPreStress<T, ConstLawPreStressType::CSL_PRE_DRIVE>  {
public:
     static constexpr ConstLawPreStressType Type = ConstLawPreStressType::CSL_PRE_DRIVE;
     static constexpr sp_grad::index_type iDim = ConstLawHelper<T>::iDim1;
     static constexpr sp_grad::index_type iNumRowsStatic = iDim;
     static_assert(ConstLawHelper<T>::iDim2 == 1);

     explicit ConstLawPreStress(const std::shared_ptr<TplDriveCaller<T>>& pTplDrv)
	  :pTplDrive(pTplDrv) {
	  DEBUGCOUTFNAME("ConstLawPreStress::ConstLawPreStress()\n");
     }

     ~ConstLawPreStress() {
	  DEBUGCOUTFNAME("ConstLawPreStress::~ConstLawPreStress()\n");
     }

     ConstLawPreStress(const ConstLawPreStress& oCslDrive)
	  :ConstLawPreStress(oCslDrive.pTplDrive) {
     }

     void Update() {
	   f0 = pTplDrive->Get();
     }

     doublereal operator()(sp_grad::index_type i) const {
	  ASSERT(i >= 1);
	  ASSERT(i <= iDim);
	  return f0(i);
     }

protected:
     T f0;
     std::shared_ptr<TplDriveCaller<T>> pTplDrive;
};

template <typename T>
class ConstLawPreStress<T, ConstLawPreStressType::CSL_PRE_CONST> {
public:
     static constexpr ConstLawPreStressType Type = ConstLawPreStressType::CSL_PRE_CONST;

     static constexpr sp_grad::index_type iDim = ConstLawHelper<T>::iDim;

     explicit ConstLawPreStress(const T& f0_a)
	  :f0(f0_a) {
     }

     ConstLawPreStress(const ConstLawPreStress& oCslDrive)
	  :f0(oCslDrive.f0) {
     }

     static void Update() {}

     doublereal operator()(sp_grad::index_type i) const {
	  ASSERT(i >= 1);
	  ASSERT(i <= iDim);
	  return f0(i);
     }

private:
     const T f0;
};

template <typename T>
class ConstLawPreStress<T, ConstLawPreStressType::CSL_PRE_NONE> {
public:
     static constexpr ConstLawPreStressType Type = ConstLawPreStressType::CSL_PRE_NONE;

     static constexpr sp_grad::index_type iDim = ConstLawHelper<T>::iDim1;

     static void Update() noexcept {}

     constexpr doublereal operator()(sp_grad::index_type) const noexcept {
	  return 0.;
     }
};

class IsotropicElasticity: public ConstitutiveLaw6D {
protected:
     IsotropicElasticity(const doublereal mu_a, const doublereal lambda_a)
	  :mu(mu_a), lambda(lambda_a) {
     }

     const doublereal mu, lambda;
};

template <typename PreStress, typename PreStrain>
class HookeanLinearElasticIsotropic: public IsotropicElasticity {
public:
     static_assert(PreStress::iDim == 6);
     static_assert(PreStrain::iDim == 6);

     HookeanLinearElasticIsotropic(const doublereal mu_a, const doublereal lambda_a, const PreStress& sigma0_a, const PreStrain& epsilon0_a)
	  :IsotropicElasticity(mu_a, lambda_a), sigma0(sigma0_a), epsilon0(epsilon0_a) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return ConstLawType::ELASTIC;
     }

     virtual ConstitutiveLaw6D* pCopy() const override {
	  HookeanLinearElasticIsotropic* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 HookeanLinearElasticIsotropic,
				 HookeanLinearElasticIsotropic(mu, lambda, sigma0, epsilon0));
	  return pCL;
     }

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
	    sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
	    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
	    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     using ConstitutiveLawAd<Vec6, Mat6x6>::Update;
     virtual void
     Update(const Vec6& Eps, const Vec6& EpsPrime) override {
	  ConstitutiveLaw6D::UpdateElasticSparse(this, Eps);
     }

     template <typename VectorType>
     void UpdateElasticTpl(const VectorType& epsilon, VectorType& sigma, const sp_grad::SpGradExpDofMapHelper<typename VectorType::ValueType>& oDofMap) {
	  using namespace sp_grad;

	  static constexpr sp_grad::index_type i1[] = {2, 1, 1};
	  static constexpr sp_grad::index_type i2[] = {3, 3, 2};

	  typedef typename VectorType::ValueType T;

	  if constexpr(std::is_same<T, doublereal>::value) {
	       sigma0.Update();
	       epsilon0.Update();
	  }

	  for (index_type i = 1; i <= 3; ++i) {
	       oDofMap.MapAssign(sigma(i), (2. * mu + lambda) * (epsilon(i) - epsilon0(i))
					 + lambda * (epsilon(i1[i - 1]) - epsilon0(i1[i - 1]) + epsilon(i2[i - 1]) - epsilon0(i2[i - 1]))
				 + sigma0(i));
	  }

	  for (index_type i = 4; i <= 6; ++i) {
	       sigma(i) = mu * (epsilon(i) - epsilon0(i)) + sigma0(i);
	  }
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  NO_OP;
     }
protected:
     PreStress sigma0;
     PreStrain epsilon0;
};

template <typename PreStress, typename PreStrain>
class HookeanLinearElasticIsotropicIncompressible: public ConstitutiveLaw7D {
public:
     using ConstitutiveLaw7D::iDimStrain;
     using ConstitutiveLaw7D::iDimStress;
     using ConstitutiveLawAd<Vec7, Mat7x7>::Update;
     static_assert(PreStress::iDim == 6);
     static_assert(PreStrain::iDim == 6);

     HookeanLinearElasticIsotropicIncompressible(const doublereal mu_a, const doublereal kappa_a, const doublereal gamma_a, const PreStress& sigma0_a, const PreStrain& epsilon0_a)
	  :mu(mu_a), kappa(kappa_a), gamma(gamma_a), sigma0(sigma0_a), epsilon0(epsilon0_a) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return ConstLawType::ELASTICINCOMPR;
     }

     virtual ConstitutiveLaw7D* pCopy() const override {
	  HookeanLinearElasticIsotropicIncompressible* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 HookeanLinearElasticIsotropicIncompressible,
				 HookeanLinearElasticIsotropicIncompressible(mu, kappa, gamma, sigma0, epsilon0));
	  return pCL;
     }

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& epsilon,
	    sp_grad::SpColVector<doublereal, iDimStress>& sigma,
	    const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap) override {
	  UpdateElasticTpl(epsilon, sigma, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& epsilon,
	    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& sigma,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap) override {
	  UpdateElasticTpl(epsilon, sigma, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& epsilon,
	    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& sigma,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap) override {
	  UpdateElasticTpl(epsilon, sigma, oDofMap);
     }

     virtual void
     Update(const Vec7& Eps, const Vec7& EpsPrime) override {
	  ConstitutiveLaw7D::UpdateElasticSparse(this, Eps);
     }

     template <typename T>
     inline void
     UpdateElasticTpl(const sp_grad::SpColVector<T, iDimStrain>& epsilon,
		      sp_grad::SpColVector<T, iDimStress>& sigma,
		      const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
	  using namespace sp_grad;

	  if constexpr(std::is_same<T, doublereal>::value) {
	       sigma0.Update();
	       epsilon0.Update();
	  }

	  const T& ptilde = epsilon(iDimStrain);

	  for (index_type i = 1; i <= 3; ++i) {
	       oDofMap.MapAssign(sigma(i), 2. * mu * (epsilon(i) - epsilon0(i)) - gamma * ptilde + sigma0(i));
	  }

	  for (index_type i = 4; i <= 6; ++i) {
	       sigma(i) = mu * (epsilon(i) - epsilon0(i)) + sigma0(i);
	  }

	  oDofMap.MapAssign(sigma(iDimStress), -(epsilon(1) + epsilon(2) + epsilon(3)
					       - (epsilon0(1) + epsilon0(2) + epsilon0(3))) - ptilde / kappa); // (pbar - ptilde) / kappa
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  NO_OP;
     }
private:
     const doublereal mu, kappa, gamma;
     PreStress sigma0;
     PreStrain epsilon0;
};

template <typename PreStress, typename PreStrain>
class HookeanLinearViscoelasticIsotropic: public HookeanLinearElasticIsotropic<PreStress, PreStrain> {
     typedef HookeanLinearElasticIsotropic<PreStress, PreStrain> BaseType;
     using BaseType::sigma0;
     using BaseType::epsilon0;
     using IsotropicElasticity::iDimStress;
     using IsotropicElasticity::iDimStrain;
     using IsotropicElasticity::lambda;
     using IsotropicElasticity::mu;
     static_assert(PreStress::iDim == 6);
     static_assert(PreStrain::iDim == 6);
public:
     HookeanLinearViscoelasticIsotropic(const doublereal mu_a, const doublereal lambda_a, const doublereal beta_a, const PreStress& sigma0_a, const PreStrain& epsilon0_a)
	  :BaseType(mu_a, lambda_a, sigma0_a, epsilon0_a), beta(beta_a) {
	  ASSERT(beta > 0.);
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return ConstLawType::VISCOELASTIC;
     }

     virtual ConstitutiveLaw6D* pCopy() const override {
	  HookeanLinearViscoelasticIsotropic* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 HookeanLinearViscoelasticIsotropic,
				 HookeanLinearViscoelasticIsotropic(mu, lambda, beta, sigma0, epsilon0));
	  return pCL;
     }

     using ConstitutiveLawAd<Vec6, Mat6x6>::Update;
     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
	    const sp_grad::SpColVector<doublereal, iDimStrain>& EpsPrime,
	    sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap) override {
	  UpdateViscoelasticTpl(Eps, EpsPrime, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
	    const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& EpsPrime,
	    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap) override {
	  UpdateViscoelasticTpl(Eps, EpsPrime, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
	    const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& EpsPrime,
	    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap) override {
	  UpdateViscoelasticTpl(Eps, EpsPrime, FTmp, oDofMap);
     }

     virtual void
     Update(const Vec6& Eps, const Vec6& EpsPrime) override {
	  ConstitutiveLaw6D::UpdateViscoelasticSparse(this, Eps, EpsPrime);
     }

     template <typename VectorType>
     void UpdateViscoelasticTpl(const VectorType& epsilon, const VectorType& epsilonP, VectorType& sigma, const sp_grad::SpGradExpDofMapHelper<typename VectorType::ValueType>& oDofMap) {
	  using namespace sp_grad;

	  static constexpr sp_grad::index_type i1[] = {2, 1, 1};
	  static constexpr sp_grad::index_type i2[] = {3, 3, 2};

	  typedef typename VectorType::ValueType T;

	  if constexpr(std::is_same<T, doublereal>::value) {
	       sigma0.Update();
	       epsilon0.Update();
	  }

	  for (index_type i = 1; i <= 3; ++i) {
	       oDofMap.MapAssign(sigma(i), (2. * mu + lambda) * ((epsilon(i) - epsilon0(i)) + beta * epsilonP(i))
				 + lambda * (epsilon(i1[i - 1]) - epsilon0(i1[i - 1]) + epsilon(i2[i - 1]) - epsilon0(i2[i - 1])
					     + beta * (epsilonP(i1[i - 1]) + epsilonP(i2[i - 1])))
				 + sigma0(i));
	  }

	  for (index_type i = 4; i <= 6; ++i) {
	       oDofMap.MapAssign(sigma(i), mu * (epsilon(i) - epsilon0(i) + beta * epsilonP(i)) + sigma0(i));
	  }
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  NO_OP;
     }
private:
     const doublereal beta;
};

class HyperElasticity {
protected:
     static void EigSym(const Mat3x3& C0, Vec3& lambda0, Mat3x3& PHI0) {
	  using namespace sp_grad;
#if defined(USE_LAPACK) && defined(HAVE_DSYEV)
	  constexpr integer N = 3;
	  constexpr integer LDA = N;
	  constexpr integer LWORK = 3 * N - 1;
	  doublereal WORK[LWORK];
	  integer INFO;

	  PHI0 = C0;

	  __FC_DECL__(dsyev)("V", "U", &N, PHI0.pGetMat(), &LDA, lambda0.pGetVec(), WORK, &LWORK, &INFO);

	  if (INFO != 0) {
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }
#else
	  C0.EigSym(lambda0, PHI0);
#endif
#ifdef DEBUG
	  {
	       const Mat3x3 LAMBDA0 = PHI0.Transpose() * C0 * PHI0;

	       const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon()) * lambda0.Norm();

	       for (index_type i = 1; i <= 3; ++i) {
		    for (index_type j = 1; j <= 3; ++j) {
			 ASSERT(fabs(LAMBDA0(i, j) - ((i == j) ? lambda0(i) : 0.)) < dTol);
		    }
	       }
	  }
#endif
     }

     template <typename VectorType, typename MatrixType, typename T>
     static void RightCauchyGreenStrainTensorNoPreStrain(const VectorType& epsilon, MatrixType& C, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
	  using namespace sp_grad;

	  static_assert(MatrixType::iNumRowsStatic == 3);
	  static_assert(MatrixType::iNumColsStatic == 3);

	  if constexpr(VectorType::iNumRowsStatic == 9) {
	       using TF = TensorFormat::Unsymmetric3x3;

	       SpMatrix<T, 3, 3> F(3, 3, 0);

	       for (index_type i = 0; i < 9; ++i) {
		    F(TF::i1[i], TF::i2[i]) = epsilon(i + 1);
	       }

	       C.MapAssign(Transpose(F) * F, oDofMap);
	  } else {
	       static_assert(VectorType::iNumRowsStatic == 6 || VectorType::iNumRowsStatic == 7);

	       static constexpr index_type i3[3][3] = {{1, 4, 6},
						       {4, 2, 5},
						       {6, 5, 3}};
	       for (index_type j = 1; j <= 3; ++j) {
		    for (index_type i = 1; i <= 3; ++i) {
			 const bool deltaij = i == j;
			 C(i, j) = (deltaij + 1.) * epsilon(i3[i - 1][j - 1]) + deltaij;
		    }
	       }
	  }
     }

     template <typename VectorType, typename ScalarType, typename PreStrain>
     static void RightCauchyGreenStrainTensorPreStrain(const VectorType& epsilon, const PreStrain& epsilon0, sp_grad::SpMatrix<ScalarType, 3, 3>& Cm, const sp_grad::SpGradExpDofMapHelper<ScalarType>& oDofMap) {
	  using namespace sp_grad;

	  Mat3x3 C0, PHI0;
	  Vec3 lambda0;

	  static_assert(PreStrain::iNumRowsStatic == 6);

	  RightCauchyGreenStrainTensorNoPreStrain(epsilon0, C0, oDofMap);

	  EigSym(C0, lambda0, PHI0);

	  for (index_type j = 1; j <= 3; ++j) {
	       if (lambda0(j) <= 0.) {
		    // physical impossible because det(F) > 0 is not fulfilled
		    throw ErrDivideByZero(MBDYN_EXCEPT_ARGS);
	       }

	       for (index_type i = 1; i <= 3; ++i) {
		    PHI0(i, j) /= sqrt(lambda0(j));
	       }
	  }

	  SpMatrixA<ScalarType, 3, 3> C;

	  RightCauchyGreenStrainTensorNoPreStrain(epsilon, C, oDofMap);

	  Cm.MapAssign(Transpose(PHI0) * C * PHI0, oDofMap);
     }

     template <typename VectorType, typename ScalarType, typename PreStrain>
     static void RightCauchyGreenStrainTensor(const VectorType& epsilon, const PreStrain& epsilon0, sp_grad::SpMatrix<ScalarType, 3, 3>& C, const sp_grad::SpGradExpDofMapHelper<ScalarType>& oDofMap) {
	  if constexpr(PreStrain::Type == ConstLawPreStressType::CSL_PRE_NONE) {
	       RightCauchyGreenStrainTensorNoPreStrain(epsilon, C, oDofMap);
	  } else {
	       RightCauchyGreenStrainTensorPreStrain(epsilon, epsilon0, C, oDofMap);
	  }
     }

     template <typename ScalarType>
     static void RightCauchyGreenStrainTensorInvariants(const sp_grad::SpMatrix<ScalarType, 3, 3>& C, sp_grad::SpMatrix<ScalarType, 3, 3>& CC, ScalarType& IC, ScalarType& IIC, ScalarType& IIIC, const sp_grad::SpGradExpDofMapHelper<ScalarType>& oDofMap) {
	  using namespace sp_grad;

	  CC.MapAssign(C * C, oDofMap);
	  oDofMap.MapAssign(IC, C(1, 1) + C(2, 2) + C(3, 3));
	  oDofMap.MapAssign(IIC, 0.5 * (IC * IC - (CC(1, 1) + CC(2, 2) + CC(3, 3))));

	  Det(C, IIIC, oDofMap);
     }
};

template <typename PreStress, typename PreStrain>
class NeoHookeanElastic: public IsotropicElasticity, public HyperElasticity  {
public:
     NeoHookeanElastic(const doublereal mu_a, const doublereal lambda_a, const PreStress& sigma0_a, const PreStrain& epsilon0_a)
	  :IsotropicElasticity(mu_a, lambda_a), sigma0(sigma0_a), epsilon0(epsilon0_a) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return ConstLawType::ELASTIC;
     }

     virtual ConstitutiveLaw6D* pCopy() const override {
	  NeoHookeanElastic* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 NeoHookeanElastic,
				 NeoHookeanElastic(mu, lambda, sigma0, epsilon0));
	  return pCL;
     }

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
	    sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
	    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
	    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     using ConstitutiveLawAd<Vec6, Mat6x6>::Update;
     virtual void
     Update(const Vec6& Eps, const Vec6& EpsPrime) override {
	  ConstitutiveLaw6D::UpdateElasticSparse(this, Eps);
     }

     template <typename VectorType>
     void UpdateElasticTpl(const VectorType& epsilon, VectorType& sigma, const sp_grad::SpGradExpDofMapHelper<typename VectorType::ValueType>& oDofMap) {
	  typedef typename VectorType::ValueType T;
	  using std::pow;
	  using namespace sp_grad;

	  // Based on Lars Kuebler 2005, chapter 2.2.1.3, page 25-26

	  SpMatrix<T, 3, 3> C(3, 3, 0), CC(3, 3, 0);
	  T IC, IIC, IIIC, gamma;

	  NeoHookeanStrainTensorAndInvariants(epsilon, C, CC, IC, IIC, IIIC, gamma, oDofMap);

	  using TFS = TensorFormat::Symmetric3x3;

	  for (index_type i = 1; i <= 6; ++i) {
	       const index_type j = TFS::i1[i - 1];
	       const index_type k = TFS::i2[i - 1];
	       const bool deltajk = (j == k);

	       oDofMap.MapAssign(sigma(i), mu * deltajk + (CC(j, k) - C(j, k) * IC + IIC * deltajk) * gamma + sigma0(i));
	  }
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  NO_OP;
     }
protected:
     template <typename T>
     void NeoHookeanStrainTensorAndInvariants(const sp_grad::SpColVector<T, 6>& epsilon, sp_grad::SpMatrix<T, 3, 3>& C, sp_grad::SpMatrix<T, 3, 3>& CC, T& IC, T& IIC, T& IIIC, T& gamma, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
	  using namespace sp_grad;

	  if constexpr(std::is_same<T, doublereal>::value) {
	      sigma0.Update();
	      epsilon0.Update();
	  }

	  RightCauchyGreenStrainTensor(epsilon, epsilon0, C, oDofMap);
	  RightCauchyGreenStrainTensorInvariants(C, CC, IC, IIC, IIIC, oDofMap);
	  oDofMap.MapAssign(gamma, (lambda * (IIIC - sqrt(IIIC)) - mu) / IIIC);
     }

     PreStress sigma0;
     PreStrain epsilon0;
};

template <typename PreStress, typename PreStrain>
class NeoHookeanViscoelastic: public NeoHookeanElastic<PreStress, PreStrain> {
     using NeoHookeanElastic<PreStress, PreStrain>::iDimStress;
     using NeoHookeanElastic<PreStress, PreStrain>::iDimStrain;
     using NeoHookeanElastic<PreStress, PreStrain>::mu;
     using NeoHookeanElastic<PreStress, PreStrain>::lambda;
     using NeoHookeanElastic<PreStress, PreStrain>::sigma0;
     using NeoHookeanElastic<PreStress, PreStrain>::epsilon0;
public:
     NeoHookeanViscoelastic(const doublereal mu_a, const doublereal lambda_a, const doublereal beta_a, const PreStress& sigma0_a, const PreStrain& epsilon0_a)
	  :NeoHookeanElastic<PreStress, PreStrain>(mu_a, lambda_a, sigma0_a, epsilon0_a), beta(beta_a) {
	  ASSERT(beta > 0.);
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return ConstLawType::VISCOELASTIC;
     }

     virtual ConstitutiveLaw6D* pCopy() const override {
	  NeoHookeanViscoelastic* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 NeoHookeanViscoelastic,
				 NeoHookeanViscoelastic(mu, lambda, beta, sigma0, epsilon0));
	  return pCL;
     }

     using ConstitutiveLawAd<Vec6, Mat6x6>::Update;
     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
	    const sp_grad::SpColVector<doublereal, iDimStress>& EpsPrime,
	    sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap) override {
	  UpdateViscoelasticTpl(Eps, EpsPrime, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
	    const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& EpsPrime,
	    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap) override {
	  UpdateViscoelasticTpl(Eps, EpsPrime, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
	    const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& EpsPrime,
	    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap) override {
	  UpdateViscoelasticTpl(Eps, EpsPrime, FTmp, oDofMap);
     }

     virtual void
     Update(const Vec6& Eps, const Vec6& EpsPrime) override {
	  ConstitutiveLaw6D::UpdateViscoelasticSparse(this, Eps, EpsPrime);
     }

     template <typename VectorType>
     void UpdateViscoelasticTpl(const VectorType& epsilon, const VectorType& epsilonP, VectorType& sigma, const sp_grad::SpGradExpDofMapHelper<typename VectorType::ValueType>& oDofMap) {
	  typedef typename VectorType::ValueType T;
	  using namespace sp_grad;
	  using TFS = TensorFormat::Symmetric3x3;

	  // Based on Lars Kuebler 2005, chapter 2.2.1.3, page 25-26
	  // In contradiction to the reference, the effect of damping
	  // is assumed proportional to initial stiffness

	  SpMatrix<T, 3, 3> C(3, 3, 0), CC(3, 3, 0);

	  T IC, IIC, IIIC, gamma;

	  this->NeoHookeanStrainTensorAndInvariants(epsilon, C, CC, IC, IIC, IIIC, gamma, oDofMap);

	  T traceGP;

	  oDofMap.MapAssign(traceGP, epsilonP(1) + epsilonP(2) + epsilonP(3));


	  for (index_type i = 1; i <= 6; ++i) {
	       const index_type j = TFS::i1[i - 1];
	       const index_type k = TFS::i2[i - 1];
	       const bool deltajk = (j == k);

	       oDofMap.MapAssign(sigma(i), mu * deltajk + (CC(j, k) - C(j, k) * IC + IIC * deltajk) * gamma
				 + (deltajk ? 2. : 1.) * mu * beta * epsilonP(i)
				 + deltajk * beta * lambda * traceGP + sigma0(i));
	  }
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  NO_OP;
     }
private:
     const doublereal beta;
};

template <typename PreStress, typename PreStrain>
class MooneyRivlinElasticBase: public HyperElasticity {
public:
     MooneyRivlinElasticBase(const doublereal C1_a, const doublereal C2_a, const doublereal kappa_a, const PreStress& sigma0_a, const PreStrain& epsilon0_a)
	  :C1(C1_a), C2(C2_a), kappa(kappa_a), sigma0(sigma0_a), epsilon0(epsilon0_a) {
     }

     template <typename T, sp_grad::index_type iDim>
     void MooneyRivlinStrainTensorAndInvariants(const sp_grad::SpColVector<T, iDim>& epsilon, sp_grad::SpMatrix<T, 3, 3>& C, sp_grad::SpMatrix<T, 3, 3>& CC, T& IC, T& IIC, T& IIIC, T& a1, T& a2, T& a3, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
	  using std::pow;
	  using namespace sp_grad;

	  if constexpr(std::is_same<T, doublereal>::value) {
	      sigma0.Update();
	      epsilon0.Update();
	  }

	  RightCauchyGreenStrainTensor(epsilon, epsilon0, C, oDofMap);
	  RightCauchyGreenStrainTensorInvariants(C, CC, IC, IIC, IIIC, oDofMap);

	  oDofMap.MapAssign(a1, 2. * C1 * pow(IIIC, -1./3.));
	  oDofMap.MapAssign(a2, 2. * C2 * pow(IIIC, -2./3.));
	  oDofMap.MapAssign(a3, a1 + a2 * IC);
     }

     template <ConstLawType::Type eConstLawType, typename T, sp_grad::index_type iDim>
     void MooneyRivlinStressTensor(const sp_grad::SpColVector<T, iDim>& epsilon, sp_grad::SpColVector<T, iDim>& sigma, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
	  using namespace sp_grad;
	  using TFS = TensorFormat::Symmetric3x3;

	  static_assert(((iDim == 6 || iDim == 9) && (eConstLawType == ConstLawType::ELASTIC)) || ((iDim == 7) && (eConstLawType == ConstLawType::ELASTICINCOMPR)));

	  // Based on K.J. Bathe
	  // Finite Element Procedures
	  // 2nd Edition
	  // Chapter 6.6.2, pages 592-594
	  // Enhancements for incompressiblity based on chapter 6.4.1, pages 561-565

	  SpMatrix<T, 3, 3> C(3, 3, 0), CC(3, 3, 0);
	  T IC, IIC, IIIC, a1, a2, a3, a4;

	  MooneyRivlinStrainTensorAndInvariants(epsilon, C, CC, IC, IIC, IIIC, a1, a2, a3, oDofMap);

	  if constexpr(eConstLawType == ConstLawType::ELASTIC) {
	       oDofMap.MapAssign(a4, (kappa * (IIIC - sqrt(IIIC)) - 2./3. * C1 * IC * pow(IIIC, -1./3.) - 4./3. * C2 * IIC * pow(IIIC, -2./3.)) / IIIC);
	  } else {
	       const T& ptilde = epsilon(iDim);
	       oDofMap.MapAssign(a4, -ptilde / sqrt(IIIC) - 2./3. * C1 * IC / pow(IIIC, 4./3.) - 4./3. * C2 * IIC / pow(IIIC, 5./3.));
	  }


	  for (index_type i = 1; i <= 6; ++i) {
	       const index_type j = TFS::i1[i - 1];
	       const index_type k = TFS::i2[i - 1];
	       const bool deltajk = (j == k);

	       oDofMap.MapAssign(sigma(i), a3 * deltajk - a2 * C(j, k) + a4 * (CC(j, k) - C(j, k) * IC + IIC * deltajk) + sigma0(i));
	  }

	  if constexpr(eConstLawType == ConstLawType::ELASTICINCOMPR) {
	       const T& ptilde = epsilon(iDim);
	       oDofMap.MapAssign(sigma(iDim), 1 - sqrt(IIIC) - ptilde / kappa); // (pbar - ptilde) / kappa
	  }

	  if constexpr(iDim == 9) {
	       using TFU = TensorFormat::Unsymmetric3x3;

	       SpMatrix<T, 3, 3> F(3, 3, 0), PK2(3, 3, 0);

	       for (index_type i = 0; i < 9; ++i) {
		    F(TFU::i1[i], TFU::i2[i]) = epsilon(i + 1);
	       }

	       for (index_type i = 0; i < 9; ++i) {
		    PK2(TFS::i1[i], TFS::i2[i]) = sigma(TFS::iv[i]);
	       }

	       const SpMatrix<T, 3, 3> PK1(F * PK2, oDofMap);

	       for (index_type i = 0; i < iDim; ++i) {
		    sigma(i + 1) = PK1(TFU::i1[i], TFU::i2[i]);
	       }
	  }
     }
protected:
     const doublereal C1, C2, kappa;
     PreStress sigma0;
     PreStrain epsilon0;
};

template <typename PreStress, typename PreStrain>
class MooneyRivlinElastic: public ConstitutiveLaw6D, private MooneyRivlinElasticBase<PreStress, PreStrain> {
     typedef MooneyRivlinElasticBase<PreStress, PreStrain> MooneyRivlinBase;
     using MooneyRivlinBase::C1;
     using MooneyRivlinBase::C2;
     using MooneyRivlinBase::kappa;
     using MooneyRivlinBase::sigma0;
     using MooneyRivlinBase::epsilon0;
public:
     MooneyRivlinElastic(const doublereal C1_a, const doublereal C2_a, const doublereal kappa_a, const PreStress& sigma0_a, const PreStrain& epsilon0_a)
	  :MooneyRivlinBase(C1_a, C2_a, kappa_a, sigma0_a, epsilon0_a) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return ConstLawType::ELASTIC;
     }

     virtual ConstitutiveLaw6D* pCopy() const override {
	  MooneyRivlinElastic* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 MooneyRivlinElastic,
				 MooneyRivlinElastic(C1, C2, kappa, sigma0, epsilon0));
	  return pCL;
     }

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
	    sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
	    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
	    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     using ConstitutiveLawAd<Vec6, Mat6x6>::Update;
     virtual void
     Update(const Vec6& Eps, const Vec6& EpsPrime) override {
	  ConstitutiveLaw6D::UpdateElasticSparse(this, Eps);
     }

     template <typename VectorType>
     void UpdateElasticTpl(const VectorType& epsilon, VectorType& sigma, const sp_grad::SpGradExpDofMapHelper<typename VectorType::ValueType>& oDofMap) {
	  this->template MooneyRivlinStressTensor<ConstLawType::ELASTIC>(epsilon, sigma, oDofMap);
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  NO_OP;
     }
};

template <typename PreStress, typename PreStrain>
class MooneyRivlinElasticDefGrad: public ConstitutiveLaw9D, private MooneyRivlinElasticBase<PreStress, PreStrain> {
     typedef MooneyRivlinElasticBase<PreStress, PreStrain> MooneyRivlinBase;
     using MooneyRivlinBase::C1;
     using MooneyRivlinBase::C2;
     using MooneyRivlinBase::kappa;
     using MooneyRivlinBase::sigma0;
     using MooneyRivlinBase::epsilon0;
public:
     MooneyRivlinElasticDefGrad(const doublereal C1_a, const doublereal C2_a, const doublereal kappa_a, const PreStress& sigma0_a, const PreStrain& epsilon0_a)
	  :MooneyRivlinBase(C1_a, C2_a, kappa_a, sigma0_a, epsilon0_a) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return ConstLawType::ELASTIC;
     }

     virtual ConstitutiveLaw9D* pCopy() const override {
	  MooneyRivlinElasticDefGrad* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 MooneyRivlinElasticDefGrad,
				 MooneyRivlinElasticDefGrad(C1, C2, kappa, sigma0, epsilon0));
	  return pCL;
     }

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
	    sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
	    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
	    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     using ConstitutiveLawAd<Vec9, Mat9x9>::Update;
     virtual void
     Update(const Vec9& Eps, const Vec9& EpsPrime) override {
	  ConstitutiveLaw9D::UpdateElasticSparse(this, Eps);
     }

     template <typename VectorType>
     void UpdateElasticTpl(const VectorType& epsilon, VectorType& sigma, const sp_grad::SpGradExpDofMapHelper<typename VectorType::ValueType>& oDofMap) {
	  this->template MooneyRivlinStressTensor<ConstLawType::ELASTIC>(epsilon, sigma, oDofMap);
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  NO_OP;
     }
};

template <typename PreStress, typename PreStrain>
class MooneyRivlinElasticIncompressible: public ConstitutiveLaw7D, private MooneyRivlinElasticBase<PreStress, PreStrain> {
     typedef MooneyRivlinElasticBase<PreStress, PreStrain> MooneyRivlinBase;
     using MooneyRivlinBase::C1;
     using MooneyRivlinBase::C2;
     using MooneyRivlinBase::kappa;
     using MooneyRivlinBase::sigma0;
     using MooneyRivlinBase::epsilon0;
public:
     using ConstitutiveLaw7D::iDimStrain;
     using ConstitutiveLaw7D::iDimStress;
     using ConstitutiveLawAd<Vec7, Mat7x7>::Update;

     MooneyRivlinElasticIncompressible(const doublereal C1_a, const doublereal C2_a, const doublereal kappa_a, const PreStress& sigma0_a, const PreStrain& epsilon0_a)
	  :MooneyRivlinBase(C1_a, C2_a, kappa_a, sigma0_a, epsilon0_a) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return ConstLawType::ELASTICINCOMPR;
     }

     virtual ConstitutiveLaw7D* pCopy() const override {
	  MooneyRivlinElasticIncompressible* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 MooneyRivlinElasticIncompressible,
				 MooneyRivlinElasticIncompressible(C1, C2, kappa, sigma0, epsilon0));
	  return pCL;
     }

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& epsilon,
	    sp_grad::SpColVector<doublereal, iDimStress>& sigma,
	    const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap) override {
	  UpdateElasticTpl(epsilon, sigma, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& epsilon,
	    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& sigma,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap) override {
	  UpdateElasticTpl(epsilon, sigma, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& epsilon,
	    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& sigma,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap) override {
	  UpdateElasticTpl(epsilon, sigma, oDofMap);
     }

     virtual void
     Update(const Vec7& Eps, const Vec7& EpsPrime) override {
	  ConstitutiveLaw7D::UpdateElasticSparse(this, Eps);
     }

     template <typename T>
     inline void
     UpdateElasticTpl(const sp_grad::SpColVector<T, iDimStrain>& epsilon,
		      sp_grad::SpColVector<T, iDimStress>& sigma,
		      const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
	  this->template MooneyRivlinStressTensor<ConstLawType::ELASTICINCOMPR>(epsilon, sigma, oDofMap);
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  NO_OP;
     }
};

struct PreStressRead {
     template <typename T>
     static std::shared_ptr<TplDriveCaller<T>> ReadPreStress(MBDynParser& HP) {
	  std::shared_ptr<TplDriveCaller<T>> sigma0(nullptr);

	  if (HP.IsKeyWord("prestress")) {
	       sigma0.reset(HP.GetTplDriveCaller<T>());
	  }

	  return sigma0;
     }

     template <typename T>
     static std::shared_ptr<TplDriveCaller<T>> ReadPreStrain(MBDynParser& HP) {
	  std::shared_ptr<TplDriveCaller<T>> epsilon0(nullptr);

	  if (HP.IsKeyWord("prestrain")) {
	       epsilon0.reset(HP.GetTplDriveCaller<T>());
	  }

	  return epsilon0;
     }

     template <typename IsotropicElasticConstLawType, typename ...Args>
     static IsotropicElasticConstLawType*
     pCreateTpl(const Args& ...args) {
	  IsotropicElasticConstLawType* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 IsotropicElasticConstLawType,
				 IsotropicElasticConstLawType(args...));

	  return pCL;
     }

     template <class ConstitutiveLawBase, template<typename, typename> class ConstitutiveLawTpl, typename ...Args>
     static ConstitutiveLawBase*
     pCreate(const std::shared_ptr<TplDriveCaller<Vec6>>& sigma0, const std::shared_ptr<TplDriveCaller<Vec6>>& epsilon0, ConstLawType::Type& CLType, const Args& ...args) {
	  ConstitutiveLawBase* pCL = nullptr;

	  typedef ConstLawPreStress<Vec6, ConstLawPreStressType::CSL_PRE_NONE> PreStressTypeNone, PreStrainTypeNone;
	  typedef ConstLawPreStress<Vec6, ConstLawPreStressType::CSL_PRE_DRIVE> PreStressTypeDrive, PreStrainTypeDrive;

	  constexpr PreStressTypeNone None;

	  // Reduce memory overhead and computational cost if pre-stress and/or pre-strain are not needed
	  if (epsilon0 && sigma0) {
	       pCL = pCreateTpl<ConstitutiveLawTpl<PreStressTypeDrive, PreStrainTypeDrive>>(args..., PreStressTypeDrive(sigma0), PreStrainTypeDrive(epsilon0));
	  } else if (sigma0) {
	       pCL = pCreateTpl<ConstitutiveLawTpl<PreStressTypeDrive, PreStrainTypeNone>>(args..., PreStressTypeDrive(sigma0), None);
	  } else if (epsilon0) {
	       pCL = pCreateTpl<ConstitutiveLawTpl<PreStressTypeNone, PreStrainTypeDrive>>(args..., None, PreStrainTypeDrive(epsilon0));
	  } else {
	       pCL = pCreateTpl<ConstitutiveLawTpl<PreStressTypeNone, PreStrainTypeNone>>(args..., None, None);
	  }

	  CLType = pCL->GetConstLawType();

	  return pCL;
     }
};

struct IsotropicElasticityRead: PreStressRead {
     struct LameParameters {
	  LameParameters(const doublereal E_, const doublereal nu_)
	       :E(E_),
		nu(nu_),
		mu(E / (2. * (1. + nu))),
		kappa(E / (3. * (1. - 2. * nu))),
		lambda(E * nu / ((1. + nu) * (1. - 2. * nu))) {
	  }

	  const doublereal E, nu, mu, kappa, lambda;
     };

     static LameParameters
     ReadLameParameters(const DataManager* pDM, MBDynParser& HP) {
	  if (!HP.IsKeyWord("E")) {
	       silent_cerr("keyword \"E\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal E = HP.GetReal(); // Young's modulus

	  if (E <= 0.) {
	       silent_cerr("E must be greater than zero at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  if (!HP.IsKeyWord("nu")) {
	       silent_cerr("keyword \"nu\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal nu = HP.GetReal(); // Poisson's ratio

	  if (nu < 0. || nu > 0.5) {
	       silent_cerr("nu must be between 0 and 0.5 at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  return LameParameters(E, nu);
     }

     template <typename T>
     static void CheckLameParameters(MBDynParser& HP, const LameParameters& lame) {
	  if constexpr(std::is_base_of<ConstitutiveLaw6D, T>::value) {
	       if (lame.nu == 0.5) {
		    silent_cerr("nu must be less than 0.5 for non incompressible 6D constitutive laws;\n"
				"use an incompressible 7D constitutive law instead at line " << HP.GetLineData() <<  "\n");
		    throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
	       }
	  }
     }

     static doublereal ReadDampingCoefficient(const DataManager* pDM, MBDynParser& HP) {
	  if (!HP.IsKeyWord("beta")) {
	       silent_cerr("keyword \"beta\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal beta = HP.GetReal(); // damping coefficient

	  if (beta < 0.) {
	       silent_cerr("beta must be greater than or equal to zero at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  return beta;
     }
};

struct HookeanLinearElasticIsotropicRead: ConstitutiveLawRead<Vec6, Mat6x6>, IsotropicElasticityRead {
     virtual ConstitutiveLaw6D*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  const auto lame = ReadLameParameters(pDM, HP);

	  CheckLameParameters<ConstitutiveLaw6D>(HP, lame);

	  const auto sigma0 = ReadPreStress<Vec6>(HP);
	  const auto epsilon0 = ReadPreStrain<Vec6>(HP);

	  return pCreate<ConstitutiveLaw6D, HookeanLinearElasticIsotropic>(sigma0, epsilon0, CLType, lame.mu, lame.lambda);
     }
};

struct HookeanLinearElasticIsotropicIncompressibleRead: ConstitutiveLawRead<Vec7, Mat7x7>, IsotropicElasticityRead {
     virtual ConstitutiveLaw7D*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  const auto lame = ReadLameParameters(pDM, HP);

	  CheckLameParameters<ConstitutiveLaw7D>(HP, lame);

	  const doublereal gamma = (3. * lame.nu) / (1. + lame.nu);

	  const auto sigma0 = ReadPreStress<Vec6>(HP);
	  const auto epsilon0 = ReadPreStrain<Vec6>(HP);

	  return pCreate<ConstitutiveLaw7D, HookeanLinearElasticIsotropicIncompressible>(sigma0, epsilon0, CLType, lame.mu, lame.kappa, gamma);
     }
};

struct HookeanLinearViscoelasticIsotropicRead: HookeanLinearElasticIsotropicRead {
     virtual ConstitutiveLaw6D*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  const auto lame = ReadLameParameters(pDM, HP);

	  CheckLameParameters<ConstitutiveLaw6D>(HP, lame);

	  const doublereal beta = ReadDampingCoefficient(pDM, HP);

	  const auto sigma0 = ReadPreStress<Vec6>(HP);
	  const auto epsilon0 = ReadPreStrain<Vec6>(HP);

	  if (!beta) {
	       return pCreate<ConstitutiveLaw6D, HookeanLinearElasticIsotropic>(sigma0, epsilon0, CLType, lame.mu, lame.lambda);
	  } else {
	       return pCreate<ConstitutiveLaw6D, HookeanLinearViscoelasticIsotropic>(sigma0, epsilon0, CLType, lame.mu, lame.lambda, beta);
	  }
     }
};

struct NeoHookeanReadElastic: ConstitutiveLawRead<Vec6, Mat6x6>, IsotropicElasticityRead {
     virtual ConstitutiveLaw6D*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  const auto lame = ReadLameParameters(pDM, HP);

	  CheckLameParameters<ConstitutiveLaw6D>(HP, lame);

	  const auto sigma0 = ReadPreStress<Vec6>(HP);
	  const auto epsilon0 = ReadPreStrain<Vec6>(HP);

	  return pCreate<ConstitutiveLaw6D, NeoHookeanElastic>(sigma0, epsilon0, CLType, lame.mu, lame.lambda);
     }
};

struct NeoHookeanReadViscoelastic: ConstitutiveLawRead<Vec6, Mat6x6>, IsotropicElasticityRead {
     virtual ConstitutiveLaw6D*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  auto lame = ReadLameParameters(pDM, HP);

	  CheckLameParameters<ConstitutiveLaw6D>(HP, lame);

	  const doublereal beta = ReadDampingCoefficient(pDM, HP);

	  const auto sigma0 = ReadPreStress<Vec6>(HP);
	  const auto epsilon0 = ReadPreStrain<Vec6>(HP);

	  if (beta == 0.) {
	       return pCreate<ConstitutiveLaw6D, NeoHookeanElastic>(sigma0, epsilon0, CLType, lame.mu, lame.lambda);
	  } else {
	       return pCreate<ConstitutiveLaw6D, NeoHookeanViscoelastic>(sigma0, epsilon0, CLType, lame.mu, lame.lambda, beta);
	  }
     }
};

struct MooneyRivlinReadElasticBase: PreStressRead {
protected:
     struct MooneyRivlinParam {
	  MooneyRivlinParam(doublereal C1_a, doublereal C2_a, doublereal kappa_a, doublereal nu_a)
	       :C1(C1_a), C2(C2_a), kappa(kappa_a), nu(nu_a) {
	  }

	  const doublereal C1, C2, kappa, nu;
     };

     static MooneyRivlinParam ReadMooneyRivlinParam(const DataManager* pDM, MBDynParser& HP) {
	  if (HP.IsKeyWord("E")) {
	       const doublereal E = HP.GetReal();

	       if (E <= 0.) {
		    silent_cerr("E must be greater than zero at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       if (!HP.IsKeyWord("nu")) {
		    silent_cerr("keyword \"nu\" expected at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const doublereal nu = HP.GetReal();

	       if (nu < 0. || nu > 0.5) {
		    silent_cerr("nu must be between zero and 0.5 at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const doublereal delta = HP.IsKeyWord("delta") ? HP.GetReal() : 0.;

	       if (delta < 0.) {
		    silent_cerr("delta must be greater than zero at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const doublereal G = E / (2. * (1. + nu));

	       const doublereal C1 = G / (2. * (1. + delta));
	       const doublereal C2 = delta * C1;
	       const doublereal kappa = E / (3. * (1. - 2. * nu));

	       return MooneyRivlinParam(C1, C2, kappa, nu);
	  } else {
	       if (!HP.IsKeyWord("C1")) {
		    silent_cerr("keyword \"C1\" expected at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const doublereal C1 = HP.GetReal();

	       if (C1 <= 0.) {
		    silent_cerr("C1 must be greater than zero at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       if (!HP.IsKeyWord("C2")) {
		    silent_cerr("keyword \"C2\" expected at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const doublereal C2 = HP.GetReal();

	       if (C2 < 0.) {
		    silent_cerr("C2 must be greater than or equal to zero at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       if (!HP.IsKeyWord("kappa")) {
		    silent_cerr("keyword \"kappa\" expected at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const doublereal kappa = HP.GetReal();

	       if (kappa <= 0.) {
		    silent_cerr("kappa must be greater than zero at line " << HP.GetLineData() << std::endl);
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const doublereal delta = C2 / C1;
	       const doublereal G = C1 * (2. * (1. + delta));
	       const doublereal nu = (3. * kappa - 2 * G) / (6 * kappa + 2 * G);

	       return MooneyRivlinParam(C1, C2, kappa, nu);
	  }
     }

     template <typename ConstitutiveLawType>
     static void CheckMooneyRivlinParam(const MooneyRivlinParam& param, MBDynParser& HP) {
	  if constexpr(!std::is_base_of<ConstitutiveLaw7D, ConstitutiveLawType>::value) {
	       if (param.nu == 0.5) {
		    silent_cerr("nu must be less than 0.5 for non incompressible 6D constitutive laws;\n"
				"An incompressible 7D constitutive law should be used instead at line "
				<< HP.GetLineData() <<  "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }
	  }
     }
};

template <template <typename, typename> class ConstitutiveLawType, typename ConstLawBaseType>
struct MooneyRivlinReadElastic: ConstitutiveLawRead<typename ConstLawBaseType::StressType, typename ConstLawBaseType::StressDerStrainType>, MooneyRivlinReadElasticBase {
     virtual ConstLawBaseType*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  const auto param = ReadMooneyRivlinParam(pDM, HP);

	  CheckMooneyRivlinParam<ConstLawBaseType>(param, HP);

	  const auto sigma0 = ReadPreStress<Vec6>(HP);
	  const auto epsilon0 = ReadPreStrain<Vec6>(HP);

	  return pCreate<ConstLawBaseType, ConstitutiveLawType>(sigma0, epsilon0, CLType, param.C1, param.C2, param.kappa);
     }
};

template <typename CSLBaseType, ConstLawType::Type CSLType, typename PreStrain>
class BilinearIsotropicHardening: public CSLBaseType {
public:
     static constexpr sp_grad::index_type iDim = CSLBaseType::iDimStress;
     static_assert(CSLBaseType::iDimStress == CSLBaseType::iDimStrain);

     static constexpr ConstLawType::Type eConstLawType = CSLType;

     BilinearIsotropicHardening(const doublereal E_a, const doublereal nu_a, const doublereal ET_a, const doublereal sigmayv_a, const PreStrain& epsilon0_a)
	  :E(E_a),
	   nu(nu_a),
	   ET(ET_a),
	   sigmayv(sigmayv_a),
	   kappa(E / (3. * (1. - 2. * nu))),
	   sigmay_prev(sigmayv),
	   sigmay_curr(sigmayv),
	   aE((1. + nu) / E),
	   EP(E * ET / (E - ET)),
	   epsilon0(epsilon0_a) {
     }

     virtual CSLBaseType* pCopy() const override {
	  BilinearIsotropicHardening* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 BilinearIsotropicHardening,
				 BilinearIsotropicHardening(E, nu, ET, sigmayv, epsilon0)); // Not considering eP_prev
	  return pCL;
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return eConstLawType;
     }

     virtual void AfterConvergence(const typename CSLBaseType::StrainType& Eps, const typename CSLBaseType::StrainType& EpsPrime) override {
	  eP_prev = eP_curr;
	  sigmay_prev = sigmay_curr;
     }

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDim>& Eps,
	    sp_grad::SpColVector<doublereal, iDim>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& Eps,
	    sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& Eps,
	    sp_grad::SpColVector<sp_grad::SpGradient, iDim>& FTmp,
	    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap) override {
	  UpdateElasticTpl(Eps, FTmp, oDofMap);
     }

     using CSLBaseType::Update;
     virtual void
     Update(const typename CSLBaseType::StrainType& Eps, const typename CSLBaseType::StrainType& EpsPrime) override {
	  CSLBaseType::UpdateElasticSparse(this, Eps);
     }

     template <typename VectorType>
     void UpdateElasticTpl(const VectorType& epsilon, VectorType& sigma, const sp_grad::SpGradExpDofMapHelper<typename VectorType::ValueType>& oDofMap) {
	  // Small strain elastoplasticity based on
	  // K.-J. Bathe Finite Element Procedures second edition 2014 ISBN 978-0-9790049-5-7
	  // Chapter 6.6.3, page 597

	  typedef typename VectorType::ValueType T;
	  using namespace sp_grad;

	  if constexpr(std::is_same<T, doublereal>::value) {
	       epsilon0.Update();
	  }

	  T em, sigmam; // mean stress, mean strain

	  oDofMap.MapAssign(em, (epsilon(1) + epsilon(2) + epsilon(3) - epsilon0(1) - epsilon0(2) - epsilon0(3)) / 3.); // equation 6.219

	  if constexpr (eConstLawType == ConstLawType::ELASTICINCOMPR) {
	       const T ptilde = epsilon(iDim);
	       sigmam = -ptilde;
	       oDofMap.MapAssign(sigma(iDim), -3. * em - ptilde / kappa); // (pbar - ptilde) / kappa
	  } else {
	       sigmam = 3. * kappa * em; // equation 6.215
	  }

	  static constexpr index_type idx1[] = {1, 2, 3, 1, 2, 3};
	  static constexpr index_type idx2[] = {1, 2, 3, 2, 3, 1};
	  static constexpr index_type idx_tens[3][3] = {{1, 4, 6},
							{4, 2, 5},
							{6, 5, 3}};

	  SpMatrix<T, 3, 3> e2(3, 3, epsilon.iGetMaxSize());

	  for (index_type i = 1; i <= 3; ++i) {
	       for (index_type j = 1; j <= 3; ++j) {
		    oDofMap.MapAssign(e2(i, j), ((i == j) ? 1. : 0.5) * (epsilon(idx_tens[i - 1][j - 1]) - epsilon0(idx_tens[i - 1][j - 1])) - em * (i == j) - eP_prev(i, j)); // equation 6.218, 6.221
	       }
	  }

	  const SpMatrix<T, 3, 3> SE(e2 / aE, oDofMap); // elastic portion of deviatoric stress tensor equation 6.239

	  T sum_e2ij_2;

	  SpGradientTraits<T>::ResizeReset(sum_e2ij_2, 0., 2 * e2.iGetMaxSize() * e2.iGetNumRows() * e2.iGetNumCols());

	  for (const auto& e2ij: e2) {
	       sum_e2ij_2 += e2ij * e2ij;
	  }

	  T d{0.};

	  if (sum_e2ij_2 > 0.) { // avoid division by zero
	       oDofMap.MapAssign(d, sqrt((3. / 2.) * sum_e2ij_2)); // equation 6.232
	  }

	  T sigma_barE; // equivalent elastic stress solution

	  oDofMap.MapAssign(sigma_barE, d / aE); // equation 6.238

	  SpMatrix<T, 3, 3> eP(3, 3, epsilon.iGetMaxSize()); // plastic strain tensor
	  SpMatrix<T, 3, 3> S(3, 3, epsilon.iGetMaxSize()); // deviatoric stress tensor

	  T sigmay{sigmay_prev}; // equivalent yield stress

	  if (sigma_barE > sigmay_prev) {
	       // yielding
	       T sigma_bar, lambda; // effective stress, yield parameter

	       oDofMap.MapAssign(sigma_bar, (2. * EP * d + 3. * sigmay_prev) / (2. * EP * aE + 3.)); // equation 6.236
	       oDofMap.MapAssign(lambda, d / sigma_bar - aE); // equation 6.230, 6.231

	       const SpMatrix<T, 3, 3> Delta_eP(SE * (aE * lambda / (aE + lambda)), oDofMap); // equation 6.240, 6.241, 6.225

	       sigmay = sigma_bar;
	       eP.MapAssign(eP_prev + Delta_eP, oDofMap);
	       S.MapAssign(SE - Delta_eP / aE, oDofMap); // stress correction equation 6.240
	  } else {
	       // not yielding
	       eP = eP_prev;
	       S = SE;
	  }

	  for (index_type i = 1; i <= 3; ++i) {
	       oDofMap.MapAssign(sigma(i), S(i, i) + sigmam); // equation 6.216
	  }

	  for (index_type i = 4; i <= 6; ++i) {
	       oDofMap.MapAssign(sigma(i), S(idx1[i - 1], idx2[i - 1])); // equation 6.216
	  }

	  UpdatePlasticStrain(eP, sigmay);
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  oData.Sync(eOwner, uLabel_a, iIndex, "eP_prev", eP_prev, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "eP_curr", eP_curr, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "sigmay_prev", sigmay_prev, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "sigmay_curr", sigmay_curr, eAction);
     }
private:
     void UpdatePlasticStrain(const sp_grad::SpMatrix<doublereal, 3, 3>& eP, doublereal sigmay) {
	  eP_curr = eP;
	  sigmay_curr = sigmay;
     }

     void UpdatePlasticStrain(const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>&,
			      const sp_grad::SpGradient&) {
     }

     void UpdatePlasticStrain(const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>&,
			      const sp_grad::GpGradProd&) {
     }

     const doublereal E, nu, ET, sigmayv, kappa;
     sp_grad::SpMatrixA<doublereal, 3, 3> eP_prev, eP_curr;
     doublereal sigmay_prev, sigmay_curr, aE, EP;
     PreStrain epsilon0;
};

template <typename PreStrain>
using BilinearIsotropicHardening6D = BilinearIsotropicHardening<ConstitutiveLaw6D, ConstLawType::ELASTIC, PreStrain>;

template <typename PreStrain>
using BilinearIsotropicHardening7D = BilinearIsotropicHardening<ConstitutiveLaw7D, ConstLawType::ELASTICINCOMPR, PreStrain>;

template <typename CSLBaseType, template<typename> class CSLType, typename CSLReadType>
struct BilinearIsotropicHardeningRead: CSLReadType, PreStressRead {
     static_assert((std::is_same<ConstitutiveLaw7D, CSLBaseType>::value && std::is_same<ConstitutiveLawRead<Vec7, Mat7x7>, CSLReadType>::value) ||
		   (std::is_same<ConstitutiveLaw6D, CSLBaseType>::value && std::is_same<ConstitutiveLawRead<Vec6, Mat6x6>, CSLReadType>::value));

     virtual CSLBaseType*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  if (!HP.IsKeyWord("E")) {
	       silent_cerr("keyword \"E\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }
	  const doublereal E = HP.GetReal(); // Young's modulus

	  if (E <= 0.) {
	       silent_cerr("E must be greater than zero at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  if (!HP.IsKeyWord("nu")) {
	       silent_cerr("keyword \"nu\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal nu = HP.GetReal(); // Poisson's ratio

	  constexpr bool bIncompressible = std::is_same<ConstitutiveLaw7D, CSLBaseType>::value;

	  if (nu < 0. || (bIncompressible ? (nu > 0.5) : (nu >= 0.5))) {
	       silent_cerr("nu must be between 0 and 0.5 at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  if (!HP.IsKeyWord("ET")) {
	       silent_cerr("keyword \"ET\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal ET = HP.GetReal(); // Plastic tangent modulus

	  if (ET < 0. || ET >= E) {
	       silent_cerr("ET must be between zero and E at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  if (!HP.IsKeyWord("sigmayv")) {
	       silent_cerr("keyword \"sigmayv\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal sigmayv = HP.GetReal(); // Equivalent initial yield stress

	  if (sigmayv <= 0.) {
	       silent_cerr("sigmayv must be greater than zero at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const auto epsilon0 = ReadPreStrain<Vec6>(HP);

	  CSLBaseType* pCL = nullptr;

	  if (epsilon0) {
	       typedef ConstLawPreStress<Vec6, ConstLawPreStressType::CSL_PRE_DRIVE> DrivePreStressType;

	       SAFENEWWITHCONSTRUCTOR(pCL,
				      CSLType<DrivePreStressType>,
				      CSLType<DrivePreStressType>(E, nu, ET, sigmayv, DrivePreStressType(epsilon0)));
	  } else {
	       typedef ConstLawPreStress<Vec6, ConstLawPreStressType::CSL_PRE_NONE> NoPreStressType;
	       constexpr NoPreStressType None;

	       SAFENEWWITHCONSTRUCTOR(pCL,
				      CSLType<NoPreStressType>,
				      CSLType<NoPreStressType>(E, nu, ET, sigmayv, None));
	  }

	  CLType = pCL->GetConstLawType();

	  return pCL;
     }
};

#if defined(USE_LAPACK) && defined(HAVE_DGETRF) && defined(HAVE_DGETRS) && defined(HAVE_DGETRI)
template <typename T, typename Tder>
class LinearViscoelasticMaxwellBase: public ConstitutiveLaw<T, Tder> {
public:
     static constexpr sp_grad::index_type iDim = ConstitutiveLaw<T, Tder>::iDimStrain;
     static_assert(ConstitutiveLaw<T, Tder>::iDimStress == ConstitutiveLaw<T, Tder>::iDimStrain); // C must be invertible

     LinearViscoelasticMaxwellBase(doublereal E0_a, const sp_grad::SpMatrix<doublereal, iDim, iDim>& C_a, DriveCaller* pTimeStepDrv_a)
	  :E0(E0_a), C(C_a), pTimeStepDrv(pTimeStepDrv_a) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
	  return ConstLawType::ELASTIC; // Because EpsPrime is not used at all
     }

protected:
     const doublereal E0;
     const sp_grad::SpMatrix<doublereal, iDim, iDim> C;
     std::unique_ptr<DriveCaller> pTimeStepDrv;
};

template <typename T, typename Tder>
class LinearViscoelasticMaxwell1: public LinearViscoelasticMaxwellBase<T, Tder> {
     typedef LinearViscoelasticMaxwellBase<T, Tder> BaseClassType;
     using BaseClassType::E0;
     using BaseClassType::C;
     using BaseClassType::pTimeStepDrv;
     using BaseClassType::F;
     using BaseClassType::FDE;
public:
     using BaseClassType::iDim;

     LinearViscoelasticMaxwell1(doublereal E0_a, doublereal E1_a, doublereal eta1_a, const sp_grad::SpMatrix<doublereal, iDim, iDim>& C_a, DriveCaller* pTimeStepDrive)
	  :BaseClassType(E0_a, C_a, pTimeStepDrive), E1(E1_a), eta1(eta1_a), EpsVPrev(iDim, 0), EpsVCurr(iDim, 0) {
	  FDE = C_a * (E0_a + E1);
     }

     virtual ConstitutiveLaw<T, Tder>* pCopy() const override {
	  LinearViscoelasticMaxwell1* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 LinearViscoelasticMaxwell1,
				 LinearViscoelasticMaxwell1(E0, E1, eta1, C, this->pTimeStepDrv->pCopy()));
	  return pCL;
     }

     using ConstitutiveLawAd<T, Tder>::Update;
     virtual void
     Update(const T& EpsCurr, const T&) override {
	  using namespace sp_grad;

	  const doublereal dt = this->pTimeStepDrv->dGet();

	  if (dt != 0.) {
	       ASSERT(dt > 0.);

	       SpMatrix<doublereal, iDim, iDim> A = C * (E1 / eta1);

	       for (index_type i = 1; i <= iDim; ++i) {
		    A(i, i) += 1. / dt;
	       }

	       EpsVCurr = C * EpsCurr * (E1 / eta1) + EpsVPrev / dt; // right hand side

	       integer INFO;
	       std::array<integer, iDim> IPIV;

	       const integer M = A.iGetNumCols();
	       const integer N = A.iGetNumRows();

	       __FC_DECL__(dgetrf)(&M, &N, &A(1, 1), &M, &IPIV[0], &INFO);

	       if (INFO != 0) {
		    silent_cerr("numeric factorization of constitutive law matrix failed with status " << INFO << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const integer NRHS = 1;

	       __FC_DECL__(dgetrs)("N", &N, &NRHS, &A(1, 1), &M, &IPIV[0], &EpsVCurr(1), &N, &INFO);

	       if (INFO != 0) {
		    silent_cerr("numeric solution of constitutive law matrix failed with status " << INFO << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       std::array<doublereal, iDim> WORK;
	       const integer LWORK = WORK.size();

	       __FC_DECL__(dgetri)(&N, &A(1, 1), &N, &IPIV[0], &WORK[0], &LWORK, &INFO);

	       if (INFO != 0) {
		    silent_cerr("numeric solution of constitutive law matrix failed with status " << INFO << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }
#ifdef DEBUG
	       const SpMatrix<doublereal, iDim, iDim> Ftmp = A * (C * (E1 / eta1) + mb_deye<Tder>(1. / dt)) - mb_deye<Tder>(1.);
	       constexpr doublereal dTolInv = 1e-10;

	       for (index_type i = 1; i <= iDim; ++i) {
		    ASSERT(Norm(Ftmp.GetCol(i)) < dTolInv);
	       }
#endif
	       const SpMatrix<doublereal, iDim, iDim> C_invA_C = C * A * C;

	       for (index_type j = 1; j <= iDim; ++j) {
		    for (index_type i = 1; i <= iDim; ++i) {
			 FDE(i, j) = C(i, j) * (E0 + E1) - C_invA_C(i, j) * (E1 * E1 / eta1);
		    }
	       }
	  } else {
	       EpsVCurr = EpsVPrev;
	       FDE = C * (E0 + E1);
	  }

	  F = SpColVector<doublereal, iDim>(C * (EpsCurr * E0 + (EpsCurr - EpsVCurr) * E1));
     }

     virtual void AfterConvergence(const T& Eps, const T&) override {
	  EpsVPrev = EpsVCurr;
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  oData.Sync(eOwner, uLabel_a, iIndex, "EpsVPrev", EpsVPrev, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "EpsVCurr", EpsVCurr, eAction);
     }
private:
     const doublereal E1, eta1;
     sp_grad::SpColVector<doublereal, iDim> EpsVPrev, EpsVCurr;
};

template <typename T, typename Tder>
class LinearViscoelasticMaxwellN: public LinearViscoelasticMaxwellBase<T, Tder> {
     typedef LinearViscoelasticMaxwellBase<T, Tder> BaseClassType;
     using BaseClassType::E0;
     using BaseClassType::C;
     using BaseClassType::pTimeStepDrv;
     using BaseClassType::F;
     using BaseClassType::FDE;
public:
     using BaseClassType::iDim;

     struct MaxwellData {
	  MaxwellData(doublereal E1_a, doublereal eta1_a)
	       :E1(E1_a), eta1(eta1_a), EpsVPrev(iDim, 0), EpsVCurr(iDim, 0) {
	  }

	  doublereal E1, eta1;
	  sp_grad::SpColVector<doublereal, iDim> EpsVPrev, EpsVCurr;
     };

     LinearViscoelasticMaxwellN(doublereal E0_a, const sp_grad::SpMatrix<doublereal, iDim, iDim>& C_a, DriveCaller* pTimeStepDrv_a, std::vector<MaxwellData>&& rgMaxwellData_a)
	  :BaseClassType(E0_a, C_a, pTimeStepDrv_a), rgMaxwellData(std::move(rgMaxwellData_a)) {
     }

     virtual ConstitutiveLaw<T, Tder>* pCopy() const override {
	  LinearViscoelasticMaxwellN* pCL = nullptr;

	  std::vector<MaxwellData> rgMaxwellDataTmp{rgMaxwellData};

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 LinearViscoelasticMaxwellN,
				 LinearViscoelasticMaxwellN(E0, C, this->pTimeStepDrv->pCopy(), std::move(rgMaxwellDataTmp)));
	  return pCL;
     }

     using ConstitutiveLawAd<T, Tder>::Update;
     virtual void
     Update(const T& EpsCurr, const T&) override {
	  using namespace sp_grad;

	  const doublereal dt = pTimeStepDrv->dGet();

	  doublereal EEq = E0;

	  for (const auto& oMaxwellData: rgMaxwellData) {
	       EEq += oMaxwellData.E1;
	  }

	  FDE = C * EEq;

	  if (dt != 0.) {
	       ASSERT(dt > 0.);

	       for (auto& oMaxwellData: rgMaxwellData) {
		    SpMatrix<doublereal, iDim, iDim> A = C * (oMaxwellData.E1 / oMaxwellData.eta1);

		    for (index_type i = 1; i <= iDim; ++i) {
			 A(i, i) += 1. / dt;
		    }

		    oMaxwellData.EpsVCurr = C * EpsCurr * (oMaxwellData.E1 / oMaxwellData.eta1) + oMaxwellData.EpsVPrev / dt; // right hand side

		    integer INFO;
		    std::array<integer, iDim> IPIV;

		    const integer M = A.iGetNumCols();
		    const integer N = A.iGetNumRows();

		    __FC_DECL__(dgetrf)(&M, &N, &A(1, 1), &M, &IPIV[0], &INFO);

		    if (INFO != 0) {
			 silent_cerr("numeric factorization of constitutive law matrix failed with status " << INFO << "\n");
			 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		    }

		    const integer NRHS = 1;

		    __FC_DECL__(dgetrs)("N", &N, &NRHS, &A(1, 1), &M, &IPIV[0], &oMaxwellData.EpsVCurr(1), &N, &INFO);

		    if (INFO != 0) {
			 silent_cerr("numeric solution of constitutive law matrix failed with status " << INFO << "\n");
			 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		    }

		    std::array<doublereal, iDim> WORK;
		    const integer LWORK = WORK.size();

		    __FC_DECL__(dgetri)(&N, &A(1, 1), &N, &IPIV[0], &WORK[0], &LWORK, &INFO);

		    if (INFO != 0) {
			 silent_cerr("numeric inversion of constitutive law matrix failed with status " << INFO << "\n");
			 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		    }
#ifdef DEBUG
		    const SpMatrix<doublereal, iDim, iDim> Ftmp = A * (C * (oMaxwellData.E1 / oMaxwellData.eta1) + mb_deye<Tder>(1. / dt)) - mb_deye<Tder>(1.);

		    constexpr doublereal dTolInv = 1e-10;

		    for (index_type i = 1; i <= iDim; ++i) {
			 ASSERT(Norm(Ftmp.GetCol(i)) < dTolInv);
		    }
#endif
		    const SpMatrix<doublereal, iDim, iDim> C_invA_C = C * A * C;
		    const doublereal d = std::pow(oMaxwellData.E1, 2)  / oMaxwellData.eta1;

		    for (index_type j = 1; j <= iDim; ++j) {
			 for (index_type i = 1; i <= iDim; ++i) {
			      FDE(i, j) -= C_invA_C(i, j) * d;
			 }
		    }
	       }
	  } else {
	       for (auto& oMaxwellData: rgMaxwellData) {
		    oMaxwellData.EpsVCurr = oMaxwellData.EpsVPrev;
	       }
	  }

	  SpColVector<doublereal, iDim> EpsEq = EpsCurr * E0;

	  for (auto& oMaxwellData: rgMaxwellData) {
	       EpsEq += (EpsCurr - oMaxwellData.EpsVCurr) * oMaxwellData.E1;
	  }

	  F = SpColVector<doublereal, iDim>(C * EpsEq);
     }

     virtual void AfterConvergence(const T& Eps, const T&) override {
	  for (auto& oMaxwellData: rgMaxwellData) {
	       oMaxwellData.EpsVPrev = oMaxwellData.EpsVCurr;
	  }
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {
	  using namespace std::string_literals;

	  for (size_t i = 0; i < rgMaxwellData.size(); ++i) {
	       const std::string strPrefix = "maxwell."s + std::to_string(i);

	       oData.Sync(eOwner, uLabel_a, iIndex, strPrefix + ".EpsVPrev", rgMaxwellData[i].EpsVPrev, eAction);
	       oData.Sync(eOwner, uLabel_a, iIndex, strPrefix + ".EpsVCurr", rgMaxwellData[i].EpsVCurr, eAction);
	  }
     }
private:
     std::vector<MaxwellData> rgMaxwellData;
};

template <typename T, typename Tder>
struct LinearViscoelasticMaxwell1Read: ConstitutiveLawRead<T, Tder> {
     virtual ConstitutiveLaw<T, Tder>*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  using namespace sp_grad;

	  if (!HP.IsKeyWord("E0")) {
	       silent_cerr("keyword \"E0\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal E0 = HP.GetReal();

	  if (E0 <= 0.) {
	       silent_cerr("E0 must be greater than zero at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  if (!HP.IsKeyWord("E1")) {
	       silent_cerr("keyword \"E1\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal E1 = HP.GetReal();

	  if (E1 < 0.) {
	       silent_cerr("E1 must be greater than or equal to zero at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  if (!HP.IsKeyWord("eta1")) {
	       silent_cerr("keyword \"eta1\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal eta1 = HP.GetReal();

	  if (eta1 <= 0.) {
	       silent_cerr("eta1 must be greater than zero at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  if (!HP.IsKeyWord("C")) {
	       silent_cerr("keyword \"C\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  typedef LinearViscoelasticMaxwell1<T, Tder> ConstLawTypeTpl;
	  constexpr index_type iDim = ConstLawTypeTpl::iDim;

	  const Tder C(HP.Get(SpMatrix<doublereal, iDim, iDim>(iDim, iDim, 0)));

	  ConstLawTypeTpl* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 ConstLawTypeTpl,
				 ConstLawTypeTpl(E0, E1, eta1, C, new TimeStepDriveCaller(pDM->pGetDrvHdl())));

	  CLType = pCL->GetConstLawType();

	  return pCL;
     }
};

template <typename T, typename Tder>
struct LinearViscoelasticMaxwellNRead: ConstitutiveLawRead<T, Tder> {
     virtual ConstitutiveLaw<T, Tder>*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  using namespace sp_grad;

	  if (!HP.IsKeyWord("E0")) {
	       silent_cerr("keyword \"E0\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const doublereal E0 = HP.GetReal();

	  if (E0 <= 0.) {
	       silent_cerr("E0 must be greater than zero at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const index_type N = HP.GetInt();

	  if (N < 1) {
	       silent_cerr("N must be greater than or equal to one at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  std::vector<typename LinearViscoelasticMaxwellN<T, Tder>::MaxwellData> rgMaxwellData;

	  rgMaxwellData.reserve(N);

	  for (index_type i = 0; i < N; ++i) {
	       if (!HP.IsKeyWord("E1")) {
		    silent_cerr("keyword \"E1\" expected at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const doublereal E1 = HP.GetReal();

	       if (E1 < 0.) {
		    silent_cerr("E1 must be greater than or equal to zero at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       if (!HP.IsKeyWord("eta1")) {
		    silent_cerr("keyword \"eta1\" expected at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       const doublereal eta1 = HP.GetReal();

	       if (eta1 <= 0.) {
		    silent_cerr("eta1 must be greater than zero at line " << HP.GetLineData() << "\n");
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       rgMaxwellData.emplace_back(E1, eta1);
	  }

	  if (!HP.IsKeyWord("C")) {
	       silent_cerr("keyword \"C\" expected at line " << HP.GetLineData() << "\n");
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  typedef LinearViscoelasticMaxwellN<T, Tder> ConstLawTypeTpl;
	  constexpr index_type iDim = ConstLawTypeTpl::iDim;

	  const Tder C(HP.Get(SpMatrix<doublereal, iDim, iDim>(iDim, iDim, 0)));

	  ConstLawTypeTpl* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 ConstLawTypeTpl,
				 ConstLawTypeTpl(E0, C, new TimeStepDriveCaller(pDM->pGetDrvHdl()), std::move(rgMaxwellData)));

	  CLType = pCL->GetConstLawType();

	  return pCL;
     }
};
#endif

#ifdef USE_MFRONT
template <typename Tstress, typename Tder, typename Tstrain = Tstress>
class MFrontGenericInterfaceCSL: public ConstitutiveLaw<Tstress, Tder, Tstrain> {
public:
     using ConstitutiveLaw<Tstress, Tder, Tstrain>::F;
     using ConstitutiveLaw<Tstress, Tder, Tstrain>::FDE;
     using ConstitutiveLaw<Tstress, Tder, Tstrain>::Update;

     explicit MFrontGenericInterfaceCSL(const mgis::behaviour::Behaviour& oBehaviourTmp, DriveCaller* pTimeStepDrv_a)
	  :oBehaviour(oBehaviourTmp),
	   oBehaviourData(oBehaviour),
	   oView(mgis::behaviour::make_view(oBehaviourData)),
	   pTimeStepDrv(pTimeStepDrv_a) {
     }

     void operator=(const MFrontGenericInterfaceCSL&)=delete;

     virtual ~MFrontGenericInterfaceCSL() {
     }

     void setMaterialProperty(const std::string_view& strName, const doublereal dValue) {
	 mgis::behaviour::setMaterialProperty(oBehaviourData.s0, strName, dValue);
	 mgis::behaviour::setMaterialProperty(oBehaviourData.s1, strName, dValue);
     }

     virtual void Update(const Tstrain& Eps, const Tstrain& EpsPrime) override {
	  mgis::behaviour::revert(oBehaviourData);

	  static constexpr size_t iNumStrain = ConstLawHelper<Tstrain>::iDim1;
	  static constexpr size_t iNumStress = ConstLawHelper<Tstress>::iDim1;

	  ASSERT(oBehaviourData.s1.gradients.size() == iNumStrain);
	  ASSERT(oBehaviourData.s1.thermodynamic_forces.size() == iNumStress);
	  ASSERT(oBehaviourData.K.size() == iNumStrain * iNumStress);
	  ASSERT(F.iGetNumRows() == iNumStress);
	  ASSERT(Eps.iGetNumRows() == iNumStrain);
	  ASSERT(FDE.iGetNumRows() == iNumStress);
	  ASSERT(FDE.iGetNumCols() == iNumStrain);

	  for (size_t i = 0; i < iNumStrain; ++i) {
	       oBehaviourData.s1.gradients[i] = dGetScale(i) * Eps(i + 1);
	  }

	  oView.dt = pTimeStepDrv->dGet();

	  // If K[0] is greater than 3.5, the consistent tangent operator must be computed.
	  oBehaviourData.K[0] = 4;

	  int status = mgis::behaviour::integrate(oView, oBehaviour);

	  if (status < 0) {
	       silent_cerr("warning: mgis::behaviour::integrate failed with status " << status << "\n");
	       throw NonlinearSolver::ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
	  }

	  for (size_t i = 0; i < iNumStress; ++i) {
	       F(i + 1) = dGetScale(i) * oBehaviourData.s1.thermodynamic_forces[i];
	  }

	  for (size_t j = 0; j < iNumStrain; ++j) {
	       for (size_t i = 0; i < iNumStress; ++i) {
		    FDE(i + 1, j + 1) = dGetScale(i) * dGetScale(j) * oBehaviourData.K[i * iNumStrain + j];
	       }
	  }

	  DEBUGLCOUT(MYDEBUG_ASSEMBLY, "Eps=" << Eps << "\n");
	  DEBUGLCOUT(MYDEBUG_ASSEMBLY, "F=" << F << "\n");
	  DEBUGLCOUT(MYDEBUG_ASSEMBLY, "FDE=" << FDE << "\n");

	  if (pedantic_out) {
	       MBDYN_LOCK_COUT();
	       mgis::behaviour::print_markdown(std::cerr, oBehaviour, oBehaviourData, 0);
	  }
     }

     virtual void AfterConvergence(const Tstrain& Eps, const Tstrain& EpsPrime) override {
	  mgis::behaviour::update(oBehaviourData);
     }

     virtual ConstLawType::Type GetConstLawType(void) const override {
	  return ConstLawType::ELASTIC;
     }

     virtual ConstitutiveLaw<Tstress, Tder, Tstrain>* pCopy() const override {
	  MFrontGenericInterfaceCSL* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 MFrontGenericInterfaceCSL,
				 MFrontGenericInterfaceCSL(oBehaviour, pTimeStepDrv->pCopy()));

	  CopyState(pCL->oBehaviourData.s0, oBehaviourData.s0);
	  CopyState(pCL->oBehaviourData.s1, oBehaviourData.s1);

	  return pCL;
     }

     virtual unsigned int iGetPrivDataIdx(const char *s) const override {
	  static constexpr struct {
	      char name[18];
	      unsigned int index;
	  } rgPrivData[] = {
	     {"requested" "time" "step", 1u},
	     {"stored" "energy", 2u},
	     {"dissipated" "energy", 3u}
	  };

	  for (const auto& oPrivData: rgPrivData) {
	      if (0 == strcmp(s, oPrivData.name)) {
		  return oPrivData.index;
	      }
	  }

	  return 0;
     }

     virtual doublereal dGetPrivData(unsigned int i) const override {
	 switch (i) {
	     case 1u:
		 return oBehaviourData.rdt;
	     case 2u:
		 return oBehaviourData.s1.stored_energy;
	     case 3:
		 return oBehaviourData.s1.dissipated_energy;
	     default:
		 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	 }
     }

     virtual unsigned int iGetNumPrivData() const override {
	  return 3u;
     }

     virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel_a, integer iIndex, RestartData::RestartAction eAction) override {

	  const size_t iNumIntStates = oBehaviourData.s0.internal_state_variables.size();
	  const size_t iNumExtStates = oBehaviourData.s0.external_state_variables.size();

	  oData.Sync(eOwner, uLabel_a, iIndex, "s0.stored_energy", oBehaviourData.s0.stored_energy, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s0.dissipated_energy", oBehaviourData.s0.dissipated_energy, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s0.mass_density", oBehaviourData.s0.mass_density, eAction);

	  oData.Sync(eOwner, uLabel_a, iIndex, "s0.gradients", oBehaviourData.s0.gradients, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s0.thermodynamic_forces", oBehaviourData.s0.thermodynamic_forces, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s0.internal_state_variables", oBehaviourData.s0.internal_state_variables, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s0.external_state_variables", oBehaviourData.s0.external_state_variables, eAction);

	  oData.Sync(eOwner, uLabel_a, iIndex, "s1.stored_energy", oBehaviourData.s1.stored_energy, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s1.dissipated_energy", oBehaviourData.s1.dissipated_energy, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s1.mass_density", oBehaviourData.s1.mass_density, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s1.gradients", oBehaviourData.s1.gradients, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s1.thermodynamic_forces", oBehaviourData.s1.thermodynamic_forces, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s1.internal_state_variables", oBehaviourData.s1.internal_state_variables, eAction);
	  oData.Sync(eOwner, uLabel_a, iIndex, "s1.external_state_variables", oBehaviourData.s1.external_state_variables, eAction);

	  static constexpr size_t iNumStrain = ConstLawHelper<Tstrain>::iDim1;
	  static constexpr size_t iNumStress = ConstLawHelper<Tstress>::iDim1;

	  if (oBehaviourData.s0.gradients.size() != iNumStrain ||
	      oBehaviourData.s0.thermodynamic_forces.size() != iNumStress ||
	      oBehaviourData.s0.internal_state_variables.size() != iNumIntStates ||
	      oBehaviourData.s0.external_state_variables.size() != iNumExtStates ||
	      oBehaviourData.s1.gradients.size() != iNumStrain ||
	      oBehaviourData.s1.thermodynamic_forces.size() != iNumStress ||
	      oBehaviourData.s1.internal_state_variables.size() != iNumIntStates ||
	      oBehaviourData.s1.external_state_variables.size() != iNumExtStates) {
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }
     }
private:
     static void CopyState(mgis::behaviour::State& oDest, const mgis::behaviour::State& oSrc) {
	  // FIXME: State::operator=(const State&) is not allowed across different behaviours!
	  oDest.material_properties = oSrc.material_properties;
	  oDest.internal_state_variables = oSrc.internal_state_variables;
	  oDest.external_state_variables = oSrc.external_state_variables;
	  oDest.gradients = oSrc.gradients;
	  oDest.thermodynamic_forces = oSrc.thermodynamic_forces;
	  oDest.dissipated_energy = oSrc.dissipated_energy;
	  oDest.mass_density = oSrc.mass_density;
     }
     static doublereal dGetScale(integer i) {
	  if constexpr(ConstLawHelper<Tstrain>::iDim1 == 6) {
	       // See: https://thelfer.github.io/tfel/web/tensors.html#symmetric-second-order-tensors
	       // s = [s11, s22, s33, sqrt(2) * s12, sqrt(2) * s13, sqrt(2) * s23]
	       constexpr doublereal c = 1. / constexpr_math::sqrt(2.);
	       return (i < 3) ? 1. : c;
	  } else {
	       static_assert(ConstLawHelper<Tstrain>::iDim1 == 9);
	       return 1.;
	  }
     }
     mgis::behaviour::Behaviour oBehaviour;
     mgis::behaviour::BehaviourData oBehaviourData;
     mgis::behaviour::BehaviourDataView oView;
     std::unique_ptr<DriveCaller> pTimeStepDrv;
};

template <typename Tstress, typename Tstrain>
struct MFrontIsSmallStrain {
     static constexpr bool value = std::is_same<Vec6, Tstress>::value && std::is_same<Vec6, Tstrain>::value;
};

template <typename Tstress, typename Tstrain>
struct MFrontIsFiniteStrain {
     static constexpr bool value = std::is_same<Vec9, Tstress>::value && std::is_same<Vec9, Tstrain>::value;
};


template <typename Tstress, typename Tstrain>
struct MFrontIsValid {
     static constexpr bool value = MFrontIsSmallStrain<Tstress, Tstrain>::value || MFrontIsFiniteStrain<Tstress, Tstrain>::value;
};

template <typename Tstress, typename Tder, typename Tstrain = Tstress, typename std::enable_if<MFrontIsValid<Tstress, Tstrain>::value, bool>::type = true>
struct MFrontGenericInterfaceCSLRead: ConstitutiveLawRead<Tstress, Tder, Tstrain> {
     virtual ConstitutiveLaw<Tstress, Tder, Tstrain>*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) override {
	  using namespace mgis;
	  using namespace mgis::behaviour;

	  static constexpr size_t iNumStrain = ConstLawHelper<Tstrain>::iDim1;

	  if (!HP.IsKeyWord("library" "path")) {
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const std::string strLibPath = HP.GetFileName();

	  if (!HP.IsKeyWord("name")) {
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }

	  const std::string strName = HP.GetStringWithDelims();

	  Behaviour oBehaviour;

	  if constexpr(iNumStrain == 6) {
	       oBehaviour = load(strLibPath, strName, Hypothesis::TRIDIMENSIONAL);
	  } else {
	       FiniteStrainBehaviourOptions oFiniteStrainOpt;
	       oFiniteStrainOpt.stress_measure = FiniteStrainBehaviourOptions::PK1;
	       oFiniteStrainOpt.tangent_operator = FiniteStrainBehaviourOptions::DPK1_DF;
	       oBehaviour = load(oFiniteStrainOpt, strLibPath, strName, Hypothesis::TRIDIMENSIONAL);
	  }

	  if (HP.IsKeyWord("parameters")) {
	       const integer iNumParam = HP.GetInt();

	       for (integer i = 0; i < iNumParam; ++i) {
		    const std::string strParamName = HP.GetStringWithDelims();
		    const doublereal dParamValue = HP.GetReal();

		    setParameter(oBehaviour, strParamName, dParamValue);
	       }
	  }

	  typedef MFrontGenericInterfaceCSL<Tstress, Tder, Tstrain> MFConstLawType;
	  MFConstLawType* pCL = nullptr;

	  SAFENEWWITHCONSTRUCTOR(pCL,
				 MFConstLawType,
				 MFConstLawType(oBehaviour, new TimeStepDriveCaller(pDM->pGetDrvHdl())));

	  if (HP.IsKeyWord("properties")) {
	      const integer iNumProp = HP.GetInt();

	      for (integer i = 0; i < iNumProp; ++i) {
		  const std::string strPropName = HP.GetStringWithDelims();
		  const doublereal dPropValue = HP.GetReal();
		  pCL->setMaterialProperty(strPropName, dPropValue);
	      }
	  }

	  CLType = pCL->GetConstLawType();

	  return pCL;
     }
};
#endif

#endif
