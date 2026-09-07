/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2025
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "automodalad.h"

AutomaticModalElemAd::AutomaticModalElemAd(const ModalNodeAd* pN)
     :AutomaticModalElem(pN),
      pModalNode(pN)
{

}

AutomaticModalElemAd::~AutomaticModalElemAd()
{
}

void AutomaticModalElemAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 6;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
AutomaticModalElemAd::AssJac(VariableSubMatrixHandler& WorkMat,
                             doublereal dCoef,
                             const VectorHandler& XCurr,
                             const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("AutomaticModalElemAd::AssJac");

     using namespace sp_grad;

     SpGradientAssVec<SpGradient>::AssJac(this,
                                          WorkMat.SetSparseGradient(),
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);

     return WorkMat;
}

SubVectorHandler&
AutomaticModalElemAd::AssRes(SubVectorHandler& WorkVec,
                             doublereal dCoef,
                             const VectorHandler& XCurr,
                             const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("AutomaticModalElemAd::AssRes");

     using namespace sp_grad;

     SpGradientAssVec<doublereal>::AssRes(this,
                                          WorkVec,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

void
AutomaticModalElemAd::AssJac(VectorHandler& JacY,
                             const VectorHandler& Y,
                             doublereal dCoef,
                             const VectorHandler& XCurr,
                             const VectorHandler& XPrimeCurr,
                             VariableSubMatrixHandler& WorkMat)
{
     DEBUGCOUTFNAME("AutomaticModalElemAd::AssJac");

     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

template <typename T>
void
AutomaticModalElemAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                             doublereal dCoef,
                             const sp_grad::SpGradientVectorHandler<T>& XCurr,
                             const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                             sp_grad::SpFunctionCall func)
{
     DEBUGCOUTFNAME("AutomaticModalElemAd::AssRes");
     using namespace sp_grad;

     const integer iRigidIndex = pModalNode->iGetFirstIndex();

     SpColVector<T, 3> xP(3, 1), g(3, 1), gP(3, 1), v(3, 1), w(3, 1);

     pModalNode->GetVCurr(xP, dCoef, func);
     pModalNode->GetgCurr(g, dCoef, func);
     pModalNode->GetgPCurr(gP, dCoef, func);
     const ::Vec3& wr = pModalNode->GetWRef();

     for (index_type i = 1; i <= 3; ++i) {
          XCurr.dGetCoef(iRigidIndex + 6 + i, v(i), dCoef);
          XCurr.dGetCoef(iRigidIndex + 9 + i, w(i), dCoef);
     }

     const SpColVector<T, 3> f1 = v - xP;
     const SpColVector<T, 3> f2 = w - MatGVec(g) * gP - MatRVec(g) * wr;

     WorkVec.AddItem(iRigidIndex + 1, f1);
     WorkVec.AddItem(iRigidIndex + 4, f2);
}
