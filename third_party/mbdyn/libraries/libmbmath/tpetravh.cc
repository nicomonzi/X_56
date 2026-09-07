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

  Tpetra port: converted from Epetra by the MBDyn project.
*/

#include "mbconfig.h"

#ifdef USE_TRILINOS
#include "matvec3.h"
#include "submat.h"
#undef HAVE_BLAS   /* avoid conflicting declaration */
#include "tpetravh.h"

/* -----------------------------------------------------------------------
 * Internal helper
 * ----------------------------------------------------------------------- */
void TpetraVectorHandler::RebuildVector(integer iNewSize, bool bZeroOut)
{
     ASSERT(iNewSize >= 0);

     const TpetraGO numGlobal = static_cast<TpetraGO>(iNewSize);
     const TpetraGO indexBase = 0;

     pMap = Teuchos::rcp(new TpetraMap(numGlobal, indexBase, pComm));
     pVec = Teuchos::rcp(new TpetraMV(pMap, bZeroOut));

     /* Sync to host and obtain a raw pointer for O(1) access. */
     //pVec->sync_host();
     oHostView = pVec->getLocalViewHost(Tpetra::Access::ReadWrite);
     pData     = oHostView.data();
}

/* -----------------------------------------------------------------------
 * Constructor / destructor
 * ----------------------------------------------------------------------- */
TpetraVectorHandler::TpetraVectorHandler(integer iSize,
                                         const Teuchos::RCP<const TpetraComm>& pComm_a)
     :pData(nullptr),
      pComm(pComm_a)
{
     RebuildVector(iSize, /*bZeroOut=*/true);
}

TpetraVectorHandler::~TpetraVectorHandler()
{
}

/* -----------------------------------------------------------------------
 * Debug
 * ----------------------------------------------------------------------- */
#ifdef DEBUG
void TpetraVectorHandler::IsValid() const
{
     ASSERT(pComm.get() != nullptr);
     ASSERT(pMap.get()  != nullptr);
     ASSERT(pVec.get()  != nullptr);
     ASSERT(pData       != nullptr);
}
#endif

/* -----------------------------------------------------------------------
 * VectorHandler interface
 * ----------------------------------------------------------------------- */
doublereal* TpetraVectorHandler::pdGetVec() const
{
#ifdef DEBUG
     IsValid();
#endif
     /* Ensure host view is up-to-date before exposing raw pointer. */
     //pVec->sync_host();
     return pData;
}

integer TpetraVectorHandler::iGetSize() const
{
#ifdef DEBUG
     IsValid();
#endif
     return static_cast<integer>(pMap->getGlobalNumElements());
}

void TpetraVectorHandler::Reset()
{
#ifdef DEBUG
     IsValid();
#endif
     pVec->putScalar(0.);
     /* Keep host view in sync. */
     //pVec->sync_host();
}

void TpetraVectorHandler::Resize(integer iNewSize)
{
#ifdef DEBUG
     IsValid();
#endif
     const integer iOldSize  = iGetSize();
     const integer iSizeCopy = std::min(iNewSize, iOldSize);

     /* Snapshot of existing data */
     std::vector<doublereal> oOld(iSizeCopy);
     //pVec->sync_host();
     for (integer i = 0; i < iSizeCopy; ++i) {
          oOld[i] = pData[i];
     }

     RebuildVector(iNewSize, /*bZeroOut=*/true);

     for (integer i = 0; i < iSizeCopy; ++i) {
          pData[i] = oOld[i];
     }
     /* Push host changes back to device. */
     //pVec->modify_host();
     //pVec->sync_device();
}

void TpetraVectorHandler::ResizeReset(integer iNewSize)
{
#ifdef DEBUG
     IsValid();
#endif
     RebuildVector(iNewSize, /*bZeroOut=*/true);
}

/* -----------------------------------------------------------------------
 * Element access – direct raw-pointer indexing (1-based MBDyn convention)
 * ----------------------------------------------------------------------- */
void TpetraVectorHandler::PutCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     pData[iRow - 1] = dCoef;
     //pVec->modify_host();
}

void TpetraVectorHandler::IncCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     pData[iRow - 1] += dCoef;
     //pVec->modify_host();
}

void TpetraVectorHandler::DecCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     pData[iRow - 1] -= dCoef;
     //pVec->modify_host();
}

const doublereal& TpetraVectorHandler::dGetCoef(integer iRow) const
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     return pData[iRow - 1];
}

const doublereal& TpetraVectorHandler::operator()(integer iRow) const
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     return pData[iRow - 1];
}

doublereal& TpetraVectorHandler::operator()(integer iRow)
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     //pVec->modify_host();
     return pData[iRow - 1];
}

/* -----------------------------------------------------------------------
 * Vec3 helpers
 * ----------------------------------------------------------------------- */
void TpetraVectorHandler::Add(integer iRow, const Vec3& v)
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow + 2 <= iGetSize());

     //pVec->modify_host();
     for (integer i = 1; i <= 3; ++i) {
          pData[iRow + i - 2] += v(i);
     }
}

void TpetraVectorHandler::Sub(integer iRow, const Vec3& v)
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow + 2 <= iGetSize());

     //pVec->modify_host();
     for (integer i = 1; i <= 3; ++i) {
          pData[iRow + i - 2] -= v(i);
     }
}

void TpetraVectorHandler::Put(integer iRow, const Vec3& v)
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow + 2 <= iGetSize());

     //pVec->modify_host();
     for (integer i = 1; i <= 3; ++i) {
          pData[iRow + i - 2] = v(i);
     }
}

/* -----------------------------------------------------------------------
 * Arithmetic operations (kept element-wise to match original semantics)
 * ----------------------------------------------------------------------- */
VectorHandler&
TpetraVectorHandler::ScalarAddMul(const VectorHandler& VH, const doublereal& d)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
#endif
     //pVec->modify_host();
     const integer iSize = iGetSize();
     for (integer i = 1; i <= iSize; ++i) {
          pData[i - 1] += d * VH(i);
     }
     return *this;
}

VectorHandler&
TpetraVectorHandler::ScalarAddMul(const VectorHandler& VH,
                                   const VectorHandler& VH1,
                                   const doublereal& d)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     VH1.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
     ASSERT(iGetSize() == VH1.iGetSize());
#endif
     //pVec->modify_host();
     const integer iSize = iGetSize();
     for (integer i = 1; i <= iSize; ++i) {
          pData[i - 1] = VH(i) + d * VH1(i);
     }
     return *this;
}

VectorHandler&
TpetraVectorHandler::ScalarMul(const VectorHandler& VH, const doublereal& d)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
#endif
     //pVec->modify_host();
     const integer iSize = iGetSize();
     for (integer i = 1; i <= iSize; ++i) {
          pData[i - 1] = d * VH(i);
     }
     return *this;
}

VectorHandler& TpetraVectorHandler::operator+=(const VectorHandler& VH)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
#endif
     //pVec->modify_host();
     const integer iSize = iGetSize();
     for (integer i = 1; i <= iSize; ++i) {
          pData[i - 1] += VH(i);
     }
     return *this;
}

VectorHandler& TpetraVectorHandler::operator+=(const SubVectorHandler& SubVH)
{
#ifdef DEBUG
     IsValid();
     SubVH.IsValid();
#endif
     //pVec->modify_host();
     SubVH.AddTo(*this);
     return *this;
}

VectorHandler& TpetraVectorHandler::operator-=(const VectorHandler& VH)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
#endif
     //pVec->modify_host();
     const integer iSize = iGetSize();
     for (integer i = 1; i <= iSize; ++i) {
          pData[i - 1] -= VH(i);
     }
     return *this;
}

VectorHandler& TpetraVectorHandler::operator*=(const doublereal& d)
{
#ifdef DEBUG
     IsValid();
#endif
     //pVec->modify_host();
     const integer iSize = iGetSize();
     for (integer i = 1; i <= iSize; ++i) {
          pData[i - 1] *= d;
     }
     return *this;
}

VectorHandler& TpetraVectorHandler::operator=(const VectorHandler& VH)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
#endif
     ResizeReset(VH.iGetSize());

     //pVec->modify_host();
     const integer iSize = iGetSize();
     for (integer i = 1; i <= iSize; ++i) {
          pData[i - 1] = VH(i);
     }
     return *this;
}

doublereal TpetraVectorHandler::Dot() const
{
#ifdef DEBUG
     IsValid();
#endif
     const integer iSize = iGetSize();
     doublereal dProd = 0.;
     for (integer i = 0; i < iSize; ++i) {
          dProd += pData[i] * pData[i];
     }
     return dProd;
}

doublereal TpetraVectorHandler::Norm() const
{
#ifdef DEBUG
     IsValid();
#endif
     return std::sqrt(Dot());
}

doublereal TpetraVectorHandler::InnerProd(const VectorHandler& VH) const
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
#endif
     const integer iSize = iGetSize();
     doublereal dProd = 0.;
     for (integer i = 1; i <= iSize; ++i) {
          dProd += pData[i - 1] * VH(i);
     }
     return dProd;
}

#endif  /* USE_TRILINOS */
