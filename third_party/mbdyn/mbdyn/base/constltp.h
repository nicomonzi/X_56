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

/* Legami costitutivi */

#ifndef CONSTLTP_H
#define CONSTLTP_H

#include "demangle.h"
#include "withlab.h"
#include "simentity.h"
#include "tpldrive.h"
#include "matvec3.h"
#include "matvec6.h"
#include "sp_gradient.h"
#include "sp_matrix_base.h"
#include "ltstrcase.h"
#include "restart_data.h"

typedef sp_grad::SpColVector<doublereal, 7> Vec7;
typedef sp_grad::SpColVector<doublereal, 9> Vec9;
typedef sp_grad::SpMatrix<doublereal, 9, 9> Mat9x9;
typedef sp_grad::SpMatrix<doublereal, 7, 7> Mat7x7;

// Required for ConstitutiveLaw7D
template <>
inline const Vec7&
mb_zero<Vec7>() {
     static const Vec7 Zero7(7, 0);

     return Zero7;
}

template <>
inline const Mat7x7&
mb_zero<Mat7x7>() {
     static const Mat7x7 Zero7x7(7, 7, 0);

     return Zero7x7;
}

template <>
inline Mat7x7
mb_deye<Mat7x7>(const doublereal d) {
     using namespace sp_grad;

     Mat7x7 A(7, 7, 0);

     for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
          A(i, i) = d;
     }

     return A;
}

template <>
inline const Vec9&
mb_zero<Vec9>() {
     static const Vec9 Zero9(9, 0);

     return Zero9;
}

// Required for ConstitutiveLaw9D
template <>
inline const Mat9x9&
mb_zero<Mat9x9>() {
     static const Mat9x9 Zero9x9(9, 9, 0);

     return Zero9x9;
}

template <>
inline Mat9x9
mb_deye<Mat9x9>(const doublereal d) = delete;

template <typename T>
struct ConstLawHelper {
     static constexpr sp_grad::index_type iDim1 = T::iNumRowsStatic;
     static constexpr sp_grad::index_type iDim2 = T::iNumColsStatic;
     static bool IsDiag(const T& FDE) { return FDE.IsDiag(); }
};

template <>
struct ConstLawHelper<doublereal> {
     static constexpr sp_grad::index_type iDim1 = 1;
     static constexpr sp_grad::index_type iDim2 = 1;
     static bool IsDiag(doublereal) { return false; } // If Eps is a scalar value, it makes no sense to call UpdateDiagonal
};

/* Tipi di cerniere deformabili */
class ConstLawType {
public:
        enum Type {
                UNKNOWN        = 0x0,

                ELASTIC        = 0x1,
                VISCOUS        = 0x2,
                INCOMPRESSIBLE = 0x4,
                VISCOELASTIC = (ELASTIC | VISCOUS),
                ELASTICINCOMPR = (ELASTIC | INCOMPRESSIBLE)
        };
};

class ConstitutiveLawCommon: public SimulationEntity, public WithLabel {
public:
        class ErrNotAvailable : public MBDynErrBase {
        public:
                ErrNotAvailable(MBDYN_EXCEPT_ARGS_DECL)
                : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU)
                {
                        silent_cerr("Constitutive law not available "
                                "for this dimensionality"
                                << std::endl);
                }
                ErrNotAvailable(std::ostream& out, MBDYN_EXCEPT_ARGS_DECL)
                : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU)
                {
                        out << "Constitutive law not available "
                                "for this dimensionality"
                                << what()
                                << std::endl;
                }
        };

        ConstitutiveLawCommon()
             :WithLabel(0) {
        }

        virtual ~ConstitutiveLawCommon() {
        }

        virtual ConstLawType::Type GetConstLawType() const = 0;

        virtual std::ostream& Restart(std::ostream& out) const {
                return out;
        }

        virtual void Restart(RestartData& oData, RestartData::RestartEntity eOwner, unsigned uLabel, integer iIndex, RestartData::RestartAction eAction);

        /* simentity */
        virtual unsigned int iGetNumDof(void) const override {
                return 0;
        }

        virtual std::ostream& DescribeDof(std::ostream& out,
                                          const char *prefix = "",
                                          bool bInitial = false) const override
        {
                return out;
        }

        virtual void DescribeDof(std::vector<std::string>& desc,
                                 bool bInitial = false, int i = -1) const override
        {
                ASSERT(i <= 0);
                desc.resize(0);
        }

        virtual std::ostream& DescribeEq(std::ostream& out,
                                         const char *prefix = "",
                                         bool bInitial = false) const override
        {
                return out;
        }

        virtual void DescribeEq(std::vector<std::string>& desc,
                                bool bInitial = false, int i = -1) const override
        {
                ASSERT(i <= 0);
                desc.resize(0);
        }

        virtual DofOrder::Order GetDofType(unsigned int i) const override {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
};

/* ConstitutiveLaw - begin */
template <typename Tstress, typename Tder, typename Tstrain = Tstress>
class ConstitutiveLawBase: public ConstitutiveLawCommon {
public:
        using SimulationEntity::AfterConvergence;
        using SimulationEntity::Update;

        typedef ErrNotAvailable Err;

        typedef Tstrain StrainType;
        typedef Tstress StressType;
        typedef Tder StressDerStrainType;
protected:
        Tstrain Epsilon;		/* strain */
        Tstrain EpsilonPrime;		/* strain rate */

        Tstress F;			/* force(stress) */
        Tder FDE;		/* force(stress) strain derivative */
        Tder FDEPrime;		/* force(stress) strain rate derivative */
public:
        static constexpr sp_grad::index_type iDimStress = ConstLawHelper<Tstress>::iDim1;
        static constexpr sp_grad::index_type iDimStrain = ConstLawHelper<Tstrain>::iDim1;

        static_assert(ConstLawHelper<Tder>::iDim1 == iDimStress);
        static_assert(ConstLawHelper<Tder>::iDim2 == iDimStrain);
        static_assert(iDimStrain >= iDimStress); // required for MFront (e.g. in order to pass the deformation gradient instead of the strain tensor)
        static_assert(ConstLawHelper<Tstress>::iDim2 == 1);
        static_assert(ConstLawHelper<Tstrain>::iDim2 == 1);

        ConstitutiveLawBase()
        : Epsilon(mb_zero<Tstrain>()), EpsilonPrime(mb_zero<Tstrain>()),
          F(mb_zero<Tstress>()), FDE(mb_zero<Tder>()), FDEPrime(mb_zero<Tder>()) {
                NO_OP;
        }

        virtual ~ConstitutiveLawBase() {
                NO_OP;
        }

        // Main interface used by all the elements without automatic differentiation
        virtual void Update(const Tstrain& Eps, const Tstrain& EpsPrime = mb_zero<Tstrain>()) = 0;

        virtual void AfterConvergence(const Tstrain& Eps, const Tstrain& EpsPrime = mb_zero<Tstrain>()) {
                NO_OP;
        }

        virtual const Tstrain& GetEpsilon(void) const {
                return Epsilon;
        }

        virtual const Tstrain& GetEpsilonPrime(void) const {
                return EpsilonPrime;
        }

        virtual const Tstress& GetF(void) const {
                return F;
        }

        virtual const Tder& GetFDE(void) const {
                return FDE;
        }

        virtual const Tder& GetFDEPrime(void) const {
                return FDEPrime;
        }
};

template <typename Tstress, typename Tder, typename Tstrain = Tstress>
class ConstitutiveLawAd: public ConstitutiveLawBase<Tstress, Tder, Tstrain> {
public:
     using ConstitutiveLawBase<Tstress, Tder, Tstrain>::iDimStress;
     using ConstitutiveLawBase<Tstress, Tder, Tstrain>::iDimStrain;
     using ConstitutiveLawBase<Tstress, Tder, Tstrain>::Update;
     using ConstitutiveLawBase<Tstress, Tder, Tstrain>::F;
     using ConstitutiveLawBase<Tstress, Tder, Tstrain>::FDE;
     using ConstitutiveLawBase<Tstress, Tder, Tstrain>::FDEPrime;

     // The following functions provide a default implementation of the interfaces,
     // used by elements with automatic differentiation, for all the constitutive laws
     // without explicit support for automatic differentiation.
     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
            sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
            sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
            sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
            const sp_grad::SpColVector<doublereal, iDimStrain>& EpsPrime,
            sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
            const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& EpsPrime,
            sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
            const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& EpsPrime,
            sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap);

protected:
     inline void
     UpdateDiagonal(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);

     inline void
     UpdateDiagonal(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
                    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
                    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap);

     inline void
     UpdateDiagonal(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                    const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& EpsPrime,
                    sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                    const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);

     inline void
     UpdateDiagonal(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
                    const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& EpsPrime,
                    sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
                    const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap);

     // Default implementation of the constitutive law matrix for constitutive laws with
     // automatic differentiation, used by elements without automatic differentiation.
     // Those functions may be instantiated by derived classes.
     template <typename ConstLaw>
     static inline void
     UpdateElasticSparse(ConstLaw* pCl, const Tstrain& Eps);

     template <typename ConstLaw>
     static inline void
     UpdateViscoelasticSparse(ConstLaw* pCl, const Tstrain& Eps, const Tstrain& EpsPrime);

     template <typename ConstLaw>
     static inline void
     UpdateElasticLocal(ConstLaw* pCLS,
                        const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                        sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                        const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);

     template <typename ConstLaw>
     static inline void
     UpdateViscoelasticLocal(ConstLaw* pCSL,
                             const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                             const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& EpsP,
                             sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                             const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);
};

template <>
class ConstitutiveLawAd<doublereal, doublereal, doublereal>: public ConstitutiveLawBase<doublereal, doublereal, doublereal> {
public:
     using ConstitutiveLawBase<doublereal, doublereal, doublereal>::iDimStress;
     using ConstitutiveLawBase<doublereal, doublereal, doublereal>::iDimStrain;
     using ConstitutiveLawBase<doublereal, doublereal, doublereal>::Update;
     using ConstitutiveLawBase<doublereal, doublereal, doublereal>::F;
     using ConstitutiveLawBase<doublereal, doublereal, doublereal>::FDE;
     using ConstitutiveLawBase<doublereal, doublereal, doublereal>::FDEPrime;

     // The following functions provide a default implementation of the interfaces,
     // used by elements with automatic differentiation, for all the constitutive laws
     // without explicit support for automatic differentiation.
     virtual void
     Update(const doublereal& Eps,
            doublereal& FTmp,
            const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap);

     virtual void
     Update(const sp_grad::SpGradient& Eps,
            sp_grad::SpGradient& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);

     virtual void
     Update(const sp_grad::GpGradProd& Eps,
            sp_grad::GpGradProd& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap);

     virtual void
     Update(const doublereal& Eps,
            const doublereal& EpsPrime,
            doublereal& FTmp,
            const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap);

     virtual void
     Update(const sp_grad::SpGradient& Eps,
            const sp_grad::SpGradient& EpsPrime,
            sp_grad::SpGradient& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);

     virtual void
     Update(const sp_grad::GpGradProd& Eps,
            const sp_grad::GpGradProd& EpsPrime,
            sp_grad::GpGradProd& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap);

protected:
     // Default implementation of the constitutive law matrix for constitutive laws with
     // automatic differentiation, used by elements without automatic differentiation.
     // Those functions may be instantiated by derived classes.
     template <typename ConstLaw>
     static void
     UpdateElasticSparse(ConstLaw* pCl, const doublereal& Eps);

     template <typename ConstLaw>
     static void
     UpdateViscoelasticSparse(ConstLaw* pCl, const doublereal& Eps, const doublereal& EpsPrime);
};

template <typename Tstress, typename Tder, typename Tstrain = Tstress>
class ConstitutiveLaw: public ConstitutiveLawAd<Tstress, Tder, Tstrain> {
     typedef ConstitutiveLawBase<Tstress, Tder, Tstrain> ConstLawBaseType;
     typedef ConstitutiveLawAd<Tstress, Tder, Tstrain> ConstLawAdType;
public:
     virtual ConstitutiveLaw* pCopy() const = 0;

     using ConstLawBaseType::iDimStress;
     using ConstLawBaseType::iDimStrain;
     using ConstLawBaseType::Update;
     using ConstLawAdType::Update;
     using ConstLawBaseType::GetConstLawType;
     using ConstLawBaseType::Restart;
     using ConstLawBaseType::AfterConvergence;
     using ConstLawBaseType::GetEpsilon;
     using ConstLawBaseType::GetEpsilonPrime;
     using ConstLawBaseType::GetF;
     using ConstLawBaseType::GetFDE;
     using ConstLawBaseType::GetFDEPrime;
     using ConstLawBaseType::iGetNumDof;
     using ConstLawBaseType::DescribeDof;
     using ConstLawBaseType::DescribeEq;
     using ConstLawBaseType::GetDofType;
     using ConstLawBaseType::F;
     using ConstLawBaseType::FDE;
     using ConstLawBaseType::FDEPrime;
};

typedef ConstitutiveLaw<doublereal, doublereal> ConstitutiveLaw1D;
typedef ConstitutiveLaw<Vec3, Mat3x3> ConstitutiveLaw3D;
typedef ConstitutiveLaw<Vec6, Mat6x6> ConstitutiveLaw6D;
typedef ConstitutiveLaw<Vec7, Mat7x7> ConstitutiveLaw7D;
typedef ConstitutiveLaw<Vec9, Mat9x9> ConstitutiveLaw9D;

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
                                                  sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
                                                  const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) == 0);

     // FIXME: Need to copy the data because vectors are not compatible
     this->Epsilon = Eps;

     this->Update(this->Epsilon, this->EpsilonPrime);

     // FIXME: Return by output argument because move constructors will not be most efficient with static vector size
     FTmp = this->GetF();
}

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                                                  sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                                                  const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) == 0);

     for (index_type i = 1; i <= iDimStress; ++i) {
          FTmp(i).ResizeReset(this->F(i), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap.GetDofMap());

          for (index_type j = 1; j <= iDimStrain; ++j) {
               Eps(j).AddDeriv(FTmp(i), this->FDE(i, j), oDofMap.GetDofMap());
          }
     }
}

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
                                                  sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
                                                  const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) == 0);

     for (index_type i = 1; i <= iDimStress; ++i) {
          FTmp(i).Reset(this->F(i));

          for (index_type j = 1; j <= iDimStrain; ++j) {
               Eps(j).InsertDeriv(FTmp(i), this->FDE(i, j));
          }
     }
}

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::Update(const sp_grad::SpColVector<doublereal, iDimStrain>& Eps,
                                                  const sp_grad::SpColVector<doublereal, iDimStrain>& EpsPrime,
                                                  sp_grad::SpColVector<doublereal, iDimStress>& FTmp,
                                                  const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) != 0);

     this->Epsilon = Eps; // FIXME: Need to copy the data because vectors are not compatible
     this->EpsilonPrime = EpsPrime;

     this->Update(this->Epsilon, this->EpsilonPrime);

     FTmp = this->GetF();
}

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                                                  const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& EpsPrime,
                                                  sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                                                  const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) != 0);

     for (index_type i = 1; i <= iDimStress; ++i) {
          FTmp(i).ResizeReset(this->F(i), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap.GetDofMap());

          for (index_type j = 1; j <= iDimStrain; ++j) {
               Eps(j).AddDeriv(FTmp(i), this->FDE(i, j), oDofMap.GetDofMap());
               EpsPrime(j).AddDeriv(FTmp(i), this->FDEPrime(i, j), oDofMap.GetDofMap());
          }
     }
}

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
                                                  const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& EpsPrime,
                                                  sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
                                                  const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) != 0);

     for (index_type i = 1; i <= iDimStress; ++i) {
          FTmp(i).Reset(this->F(i));

          for (index_type j = 1; j <= iDimStrain; ++j) {
               Eps(j).InsertDeriv(FTmp(i), this->FDE(i, j));
               EpsPrime(j).InsertDeriv(FTmp(i), this->FDEPrime(i, j));
          }
     }
}

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::UpdateDiagonal(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                                                          sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                                                          const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) == 0);
     ASSERT(ConstLawHelper<Tder>::IsDiag(this->FDE));

     for (index_type i = 1; i <= iDimStress; ++i) {
          FTmp(i).ResizeReset(this->F(i), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap.GetDofMap());
          Eps(i).AddDeriv(FTmp(i), this->FDE(i, i), oDofMap.GetDofMap());
     }
}

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::UpdateDiagonal(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
                                                          sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
                                                          const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) == 0);
     ASSERT(ConstLawHelper<Tder>::IsDiag(this->FDE));

     for (index_type i = 1; i <= iDimStress; ++i) {
          FTmp(i).Reset(this->F(i));
          Eps(i).InsertDeriv(FTmp(i), this->FDE(i, i));
     }
}

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::UpdateDiagonal(const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                                                          const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& EpsPrime,
                                                          sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                                                          const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) != 0);
     ASSERT(ConstLawHelper<Tder>::IsDiag(this->FDE));
     ASSERT(ConstLawHelper<Tder>::IsDiag(this->FDEPrime));

     for (index_type i = 1; i <= iDimStress; ++i) {
          FTmp(i).ResizeReset(this->F(i), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap.GetDofMap());
          Eps(i).AddDeriv(FTmp(i), this->FDE(i, i), oDofMap.GetDofMap());
          EpsPrime(i).AddDeriv(FTmp(i), this->FDEPrime(i, i), oDofMap.GetDofMap());
     }
}

template <typename Tstress, typename Tder, typename Tstrain>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::UpdateDiagonal(const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& Eps,
                                                          const sp_grad::SpColVector<sp_grad::GpGradProd, iDimStrain>& EpsPrime,
                                                          sp_grad::SpColVector<sp_grad::GpGradProd, iDimStress>& FTmp,
                                                          const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) != 0);
     ASSERT(ConstLawHelper<Tder>::IsDiag(this->FDE));
     ASSERT(ConstLawHelper<Tder>::IsDiag(this->FDEPrime));

     for (index_type i = 1; i <= iDimStress; ++i) {
          FTmp(i).Reset(this->F(i));
          Eps(i).InsertDeriv(FTmp(i), this->FDE(i, i));
          EpsPrime(i).InsertDeriv(FTmp(i), this->FDEPrime(i, i));
     }
}

template <typename Tstress, typename Tder, typename Tstrain>
template <typename ConstLaw>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::UpdateViscoelasticSparse(ConstLaw* pCl, const Tstrain& Eps, const Tstrain& EpsPrime)
{
        using namespace sp_grad;
        SpColVectorA<SpGradient, iDimStrain> gEps, gEpsPrime;
        SpColVectorA<SpGradient, iDimStress> gF;

        pCl->Epsilon = Eps;
        pCl->EpsilonPrime = EpsPrime;

        for (index_type i = 1; i <= iDimStrain; ++i)
        {
             gEps(i).Reset(Eps(i), i, 1.);
             gEpsPrime(i).Reset(EpsPrime(i), i + iDimStrain, 1.);
        }

        SpGradExpDofMapHelper<SpGradient> oDofMap;

        oDofMap.GetDofStat(gEps);
        oDofMap.GetDofStat(gEpsPrime);
        oDofMap.Reset();
        oDofMap.InsertDof(gEps);
        oDofMap.InsertDof(gEpsPrime);
        oDofMap.InsertDone();

        pCl->UpdateViscoelasticTpl(gEps, gEpsPrime, gF, oDofMap);

        for (index_type j = 1; j <= iDimStrain; ++j) {
             for (index_type i = 1; i <= iDimStress; ++i) {
                  pCl->FDE(i, j) = pCl->FDEPrime(i, j) = 0.;
             }
        }

        for (index_type i = 1; i <= iDimStress; ++i) {
                pCl->F(i) = gF(i).dGetValue();

                for (const auto& oDer: gF(i)) {
                     ASSERT(oDer.iDof >= 1);
                     ASSERT(oDer.iDof <= 2 * iDimStrain);

                     if (oDer.iDof <= iDimStrain) {
                          pCl->FDE(i, oDer.iDof) += oDer.dDer;
                     } else {
                          pCl->FDEPrime(i, oDer.iDof - iDimStrain) += oDer.dDer;
                     }
                }
        }
}

template <typename T, typename Tder, typename Tstrain>
template <typename ConstLaw>
void
ConstitutiveLawAd<T, Tder, Tstrain>::UpdateElasticSparse(ConstLaw* pCl, const Tstrain& Eps)
{
        using namespace sp_grad;
        SpColVectorA<SpGradient, iDimStrain> gEps;
        SpColVectorA<SpGradient, iDimStress> gF;

        pCl->Epsilon = Eps;

        for (index_type i = 1; i <= iDimStrain; ++i)
        {
             gEps(i).Reset(Eps(i), i, 1.);
        }

        SpGradExpDofMapHelper<SpGradient> oDofMap;

        oDofMap.GetDofStat(gEps);
        oDofMap.Reset();
        oDofMap.InsertDof(gEps);
        oDofMap.InsertDone();

        pCl->UpdateElasticTpl(gEps, gF, oDofMap);

        for (index_type j = 1; j <= iDimStrain; ++j) {
             for (index_type i = 1; i <= iDimStress; ++i) {
                  pCl->FDE(i, j) = 0.;
             }
        }

        for (index_type i = 1; i <= iDimStress; ++i) {
                pCl->F(i) = gF(i).dGetValue();

                for (const auto& oDer: gF(i)) {
                     ASSERT(oDer.iDof >= 1);
                     ASSERT(oDer.iDof <= iDimStrain);
                     pCl->FDE(i, oDer.iDof) += oDer.dDer;
                }
        }
}

template <typename Tstress, typename Tder, typename Tstrain>
template <typename ConstLaw>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::UpdateElasticLocal(ConstLaw* pCSL,
                                                              const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                                                              sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                                                              const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap)
{
     using namespace sp_grad;

     constexpr SpGradExpDofMapHelper<GpGradProd> oDofMapLocal;

     SpColVector<GpGradProd, iDimStrain> gEps(iDimStrain, 0);
     SpColVector<GpGradProd, iDimStress> gF(iDimStress, 0);

     for (index_type j = 1; j <= iDimStrain; ++j) {
          for (index_type i = 1; i <= iDimStrain; ++i) {
               const doublereal delta_ij = i == j ? 1. : 0.;

               gEps(i).Reset(Eps(i).dGetValue(), delta_ij);
          }

          pCSL->UpdateElasticTpl(gEps, gF, oDofMapLocal);

          for (index_type i = 1; i <= iDimStress; ++i) {
               pCSL->FDE(i, j) = gF(i).dGetDer();
          }
     }

     ASSERT(pCSL->iGetNumDof() == 0);
     ASSERT((pCSL->GetConstLawType() & ConstLawType::VISCOUS) == 0);

     for (index_type i = 1; i <= iDimStress; ++i) {
          ASSERT(std::fabs(gF(i).dGetValue() - pCSL->F(i)) < sqrt(std::numeric_limits<doublereal>::epsilon()) * std::max(1., std::fabs(pCSL->F(i))));

          FTmp(i).ResizeReset(gF(i).dGetValue(), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap.GetDofMap());

          for (index_type j = 1; j <= iDimStrain; ++j) {
               Eps(j).AddDeriv(FTmp(i), pCSL->FDE(i, j), oDofMap.GetDofMap());
          }
     }
}

template <typename Tstress, typename Tder, typename Tstrain>
template <typename ConstLaw>
void
ConstitutiveLawAd<Tstress, Tder, Tstrain>::UpdateViscoelasticLocal(ConstLaw* pCSL,
                                                                   const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& Eps,
                                                                   const sp_grad::SpColVector<sp_grad::SpGradient, iDimStrain>& EpsP,
                                                                   sp_grad::SpColVector<sp_grad::SpGradient, iDimStress>& FTmp,
                                                                   const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap)
{
     using namespace sp_grad;

     constexpr SpGradExpDofMapHelper<GpGradProd> oDofMapLocal;

     SpColVector<GpGradProd, iDimStrain> gEps(iDimStrain, 0), gEpsP(iDimStrain, 0);
     SpColVector<GpGradProd, iDimStress> gF(iDimStress, 0);

     for (index_type j = 1; j <= iDimStrain; ++j) {
          for (index_type i = 1; i <= iDimStrain; ++i) {
               const doublereal delta_ij = i == j ? 1. : 0.;

               gEps(i).Reset(Eps(i).dGetValue(), delta_ij);
               gEpsP(i).Reset(EpsP(i).dGetValue(), 0.);
          }

          pCSL->UpdateViscoelasticTpl(gEps, gEpsP, gF, oDofMapLocal);

          for (index_type i = 1; i <= iDimStress; ++i) {
               pCSL->FDE(i, j) = gF(i).dGetDer();
          }
     }

     for (index_type j = 1; j <= iDimStrain; ++j) {
          for (index_type i = 1; i <= iDimStrain; ++i) {
               const doublereal delta_ij = i == j ? 1. : 0.;

               gEps(i).Reset(Eps(i).dGetValue(), 0.);
               gEpsP(i).Reset(EpsP(i).dGetValue(), delta_ij);
          }

          pCSL->UpdateViscoelasticTpl(gEps, gEpsP, gF, oDofMapLocal);

          for (index_type i = 1; i <= iDimStress; ++i) {
               pCSL->FDEPrime(i, j) = gF(i).dGetDer();
          }
     }

     ASSERT(pCSL->iGetNumDof() == 0);
     ASSERT((pCSL->GetConstLawType() & ConstLawType::VISCOUS) == 0);

     for (index_type i = 1; i <= iDimStress; ++i) {
          ASSERT(std::fabs(gF(i).dGetValue() - pCSL->F(i)) < sqrt(std::numeric_limits<doublereal>::epsilon()) * std::max(1., std::fabs(pCSL->F(i))));

          FTmp(i).ResizeReset(gF(i).dGetValue(), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap.GetDofMap());

          for (index_type j = 1; j <= iDimStrain; ++j) {
               Eps(j).AddDeriv(FTmp(i), pCSL->FDE(i, j), oDofMap.GetDofMap());
               EpsP(j).AddDeriv(FTmp(i), pCSL->FDEPrime(i, j), oDofMap.GetDofMap());
          }
     }
}

template <typename ConstLaw>
void
ConstitutiveLawAd<doublereal, doublereal, doublereal>::UpdateElasticSparse(ConstLaw* pCl, const doublereal& Eps)
{
        using namespace sp_grad;

        constexpr SpGradExpDofMapHelper<GpGradProd> oDofMap;

        GpGradProd gEps(Eps, 1.);
        GpGradProd gF;

        pCl->UpdateElasticTpl(gEps, gF, oDofMap);

        pCl->Epsilon = Eps;

        pCl->F = gF.dGetValue();

        pCl->FDE = gF.dGetDeriv();
}

/* ConstitutiveLaw - end */


/* ConstitutiveLawOwner - begin */

// template <class T, class Tder>
// class ConstitutiveLawOwner : public SimulationEntity {
// protected:
//         mutable ConstitutiveLaw<T, Tder>* pConstLaw;
// 
// public:
//         ConstitutiveLawOwner(const ConstitutiveLaw<T, Tder>* pCL)
//         : pConstLaw(const_cast<ConstitutiveLaw<T, Tder> *>(pCL)) {
//                 ASSERT(pCL != NULL);
//         };
// 
//         virtual ~ConstitutiveLawOwner(void) {
//                 ASSERT(pConstLaw != NULL);
//                 if (pConstLaw != NULL) {
//                         SAFEDELETE(pConstLaw);
//                 }
//         };
// 
//         inline ConstitutiveLaw<T, Tder>* pGetConstLaw(void) const {
//                 ASSERT(pConstLaw != NULL);
//                 return pConstLaw;
//         };
// 
// 	using SimulationEntity::Update;
//         inline void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
//                 ASSERT(pConstLaw != NULL);
//                 pConstLaw->Update(Eps, EpsPrime);
//         };
// 
// 	using SimulationEntity::AfterConvergence;
//         inline void AfterConvergence(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
//                 ASSERT(pConstLaw != NULL);
//                 pConstLaw->AfterConvergence(Eps, EpsPrime);
//         };
// 
//         inline const T& GetF(void) const {
//                 ASSERT(pConstLaw != NULL);
//                 return pConstLaw->GetF();
//         };
// 
//         inline const Tder& GetFDE(void) const {
//                 ASSERT(pConstLaw != NULL);
//                 return pConstLaw->GetFDE();
//         };
// 
//         inline const Tder& GetFDEPrime(void) const {
//                 ASSERT(pConstLaw != NULL);
//                 return pConstLaw->GetFDEPrime();
//         };
// 
//         /* simentity */
//         virtual unsigned int iGetNumDof(void) const override {
//                 ASSERT(pConstLaw != NULL);
//                 return pConstLaw->iGetNumDof();
//         };
// 
//         virtual std::ostream& DescribeDof(std::ostream& out,
//                         const char *prefix = "",
//                         bool bInitial = false) const override
//         {
//                 return out;
//         };
// 
//         virtual void DescribeDof(std::vector<std::string>& desc,
//                         bool bInitial = false, int i = -1) const override
//         {
//                 ASSERT(i <= 0);
//                 desc.resize(0);
//         };
// 
//         virtual std::ostream& DescribeEq(std::ostream& out,
//                         const char *prefix = "",
//                         bool bInitial = false) const override
//         {
//                 return out;
//         };
// 
//         virtual void DescribeEq(std::vector<std::string>& desc,
//                         bool bInitial = false, int i = -1) const override
//         {
//                 ASSERT(i <= 0);
//                 desc.resize(0);
//         };
// 
//         virtual DofOrder::Order GetDofType(unsigned int i) const override {
//                 ASSERT(pConstLaw != NULL);
//                 return pConstLaw->GetDofType(i);
//         };
// 
//         /*
//          * Metodi per l'estrazione di dati "privati".
//          * Si suppone che l'estrattore li sappia interpretare.
//          * Come default non ci sono dati privati estraibili
//          */
//         virtual unsigned int iGetNumPrivData(void) const override {
//                 return pConstLaw->iGetNumPrivData();
//         };
// 
//         /*
//          * Maps a string (possibly with substrings) to a private data;
//          * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0
//          * in case of unrecognized data; error must be handled by caller
//          */
//         virtual unsigned int iGetPrivDataIdx(const char *s) const override {
//                 return pConstLaw->iGetPrivDataIdx(s);
//         };
// 
//         /*
//          * Returns the current value of a private data
//          * with 0 < i <= iGetNumPrivData()
//          */
//         virtual doublereal dGetPrivData(unsigned int i) const override {
//                 return pConstLaw->dGetPrivData(i);
//         };
// 
//         virtual std::ostream& OutputAppend(std::ostream& out) const override {
//                 return pConstLaw->OutputAppend(out);
//         };
// 
//         virtual void NetCDFOutputAppend(OutputHandler& OH) const override {
//                 return pConstLaw->NetCDFOutputAppend(OH);
//         };
// 
//         virtual void OutputAppendPrepare(OutputHandler& OH, const std::string& name) override {
//                 pConstLaw->OutputAppendPrepare(OH, name);
//         };
// 
// };
// 
// typedef ConstitutiveLawOwner<doublereal, doublereal> ConstitutiveLaw1DOwner;
// typedef ConstitutiveLawOwner<Vec3, Mat3x3> ConstitutiveLaw3DOwner;
// typedef ConstitutiveLawOwner<Vec6, Mat6x6> ConstitutiveLaw6DOwner;

/* ConstitutiveLawOwner - end */

/* prototype of the template functional object: reads a constitutive law */
template <typename Tstress, typename Tder, typename Tstrain = Tstress>
struct ConstitutiveLawRead {
        virtual ~ConstitutiveLawRead( void ) { NO_OP; };
        virtual ConstitutiveLaw<Tstress, Tder, Tstrain>*
        Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) = 0;
};

/* functions that read a constitutive law */
ConstitutiveLaw<doublereal, doublereal> *
ReadCL1D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);

ConstitutiveLaw<Vec3, Mat3x3> *
ReadCL3D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);

ConstitutiveLaw<Vec6, Mat6x6> *
ReadCL6D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);

ConstitutiveLaw<Vec7, Mat7x7> *
ReadCL7D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);

ConstitutiveLaw<Vec9, Mat9x9>*
ReadCL9D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);

/* constitutive law registration functions: call to register one */
bool
SetCL1D(const std::string& name, ConstitutiveLawRead<doublereal, doublereal>* rf);

bool
SetCL3D(const std::string& name, ConstitutiveLawRead<Vec3, Mat3x3>* rf);

bool
SetCL6D(const std::string& name, ConstitutiveLawRead<Vec6, Mat6x6>* rf);

bool
SetCL7D(const std::string& name, ConstitutiveLawRead<Vec7, Mat7x7>* rf);

bool
SetCL9D(const std::string& name, ConstitutiveLawRead<Vec9, Mat9x9>* rf);

void ClearCL();

#endif /* CONSTLTP_H */
