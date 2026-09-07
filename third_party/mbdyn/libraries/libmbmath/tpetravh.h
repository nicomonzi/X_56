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

#ifndef ___TPETRA_VECTOR_HANDLER__INCLUDED___
#define ___TPETRA_VECTOR_HANDLER__INCLUDED___

#ifdef USE_TRILINOS
#include "vh.h"
#include "tpetra_types.h"

class TpetraVectorHandler: public VectorHandler {
public:
     TpetraVectorHandler(integer iSize,
                         const Teuchos::RCP<const TpetraComm>& pComm);
     virtual ~TpetraVectorHandler();

#ifdef DEBUG
     virtual void IsValid() const;
#endif
     virtual doublereal* pdGetVec() const;

     virtual integer iGetSize() const;

     virtual void Reset();

     virtual void Resize(integer iNewSize);

     virtual void ResizeReset(integer iNewSize);

     virtual void PutCoef(integer iRow, const doublereal& dCoef);

     virtual void IncCoef(integer iRow, const doublereal& dCoef);

     virtual void DecCoef(integer iRow, const doublereal& dCoef);

     virtual const doublereal& dGetCoef(integer iRow) const;

     virtual const doublereal& operator()(integer iRow) const;

     virtual doublereal& operator()(integer iRow);

     virtual void Add(integer iRow, const Vec3& v);

     virtual void Sub(integer iRow, const Vec3& v);

     virtual void Put(integer iRow, const Vec3& v);

     virtual VectorHandler&
     ScalarAddMul(const VectorHandler& VH, const doublereal& d);

     virtual VectorHandler&
     ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
                  const doublereal& d);

     virtual VectorHandler&
     ScalarMul(const VectorHandler& VH, const doublereal& d);

     virtual VectorHandler& operator+=(const VectorHandler& VH);

     virtual VectorHandler& operator+=(const SubVectorHandler& SubVH);

     virtual VectorHandler& operator-=(const VectorHandler& VH);

     virtual VectorHandler& operator*=(const doublereal& d);

     virtual VectorHandler& operator=(const VectorHandler& VH);

     virtual doublereal Dot() const;

     virtual doublereal Norm() const;

     virtual doublereal InnerProd(const VectorHandler& VH) const;

     Teuchos::RCP<const TpetraMV> pGetTpetraVector() const { return pVec; }
     Teuchos::RCP<TpetraMV>       pGetTpetraVector()       { return pVec; }

private:
     /* Rebuild the map and vector with a new size, preserving the communicator. */
     void RebuildVector(integer iNewSize, bool bZeroOut);

     /*
      * Host-side 1-D view into the Tpetra vector data.
      * We keep a raw pointer for O(1) element access so that the
      * per-element VectorHandler methods stay as cheap as before.
      * The pointer is refreshed whenever the vector is rebuilt.
      */
     doublereal* pData;

     Teuchos::RCP<const TpetraComm>   pComm;
     Teuchos::RCP<const TpetraMap>    pMap;
     Teuchos::RCP<TpetraMV>       pVec;

     /* Kokkos host view – kept alive so the pointer stays valid. */
     mutable TpetraMV::dual_view_type::t_host oHostView;
};

#endif  /* USE_TRILINOS */
#endif  /* ___TPETRA_VECTOR_HANDLER__INCLUDED___ */
