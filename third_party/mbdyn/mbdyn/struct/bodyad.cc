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
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2023) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include "bodyad.h"

template<class B, class N>
BodyAd<B, N>::BodyAd(unsigned int uL, const N *pNode_a,
            doublereal dMassTmp, const Vec3& XgcTmp, const Mat3x3& JTmp_a,
            flag fOut)
     : B(uL, pNode_a, dMassTmp, XgcTmp, JTmp_a, fOut),
     pNode(pNode_a)
{
     ASSERT(pNode != nullptr);
     ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
     ASSERT(dMassTmp > 0.);
}

template<class B, class N>
BodyAd<B, N>::~BodyAd(void)
{
}

template <class B, class N> 
template<typename T>
void
BodyAd<B, N>::AssVecRBK_int(const RigidBodyKinematics* const pRBK,
                      const sp_grad::SpColVector<T, 3>& X,
                      const sp_grad::SpColVector<T, 3>& V,
                      const sp_grad::SpColVector<T, 3>& W,
                      const sp_grad::SpColVector<T, 3>& STmp_a,
                      const sp_grad::SpMatrix<T, 3, 3>& JTmp_a,
                      sp_grad::SpGradientAssVec<T>& WorkVec,
                      const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     const SpColVector<T, 3> s0 = X * B::dMass + STmp_a;

     // force
     const SpColVector<T, 3> F(pRBK->GetXPP() * -B::dMass
                               - Cross(pRBK->GetWP(), s0, oDofMap)
                               - Cross(pRBK->GetW(), Cross(pRBK->GetW(), s0, oDofMap), oDofMap), oDofMap);

     const SpColVector<T, 3> a(pRBK->GetXPP()
                               + Cross(pRBK->GetWP(), X, oDofMap)
                               + Cross(pRBK->GetW(), Cross(pRBK->GetW(), X, oDofMap), oDofMap)
                               + Cross(pRBK->GetW(), V, oDofMap), oDofMap);

     const SpColVector<T, 3> JTmpWRBK(JTmp_a * pRBK->GetW(), oDofMap);

     // moment
     const SpColVector<T, 3> M(-Cross(STmp_a, a, oDofMap)
                               - Cross(pRBK->GetW(), JTmpWRBK, oDofMap)
                               - SpColVector<T, 3>(JTmp_a * pRBK->GetWP(), oDofMap)
                               - Cross(W, JTmpWRBK, oDofMap)
                               + SpColVector<T, 3>(JTmp_a * Cross(W, pRBK->GetW(), oDofMap), oDofMap)
                               - Cross(V, Cross(pRBK->GetW(), STmp_a, oDofMap), oDofMap), oDofMap);

     const integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();

     WorkVec.AddItem(iFirstMomentumIndex + 1, F);
     WorkVec.AddItem(iFirstMomentumIndex + 4, M);
}

template<class B, class N>
void BodyAd<B, N>::UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& S,
                           const sp_grad::SpMatrix<doublereal, 3, 3>& J) const
{
     for (integer i = 1; i <= 3; ++i) {
          B::STmp(i) = S(i);
     }

     for (integer j = 1; j <= 3; ++j) {
          for (integer i = 1; i <= 3; ++i) {
               B::JTmp(i, j) = J(i, j);
          }
     }
}

DynamicBodyAd::DynamicBodyAd(unsigned int uL,
                             const DynamicStructNodeAd* pN,
                             doublereal dMass_a,
                             const Vec3& Xgc_a,
                             const Mat3x3& J,
                             flag fOut)
     : BodyAd<DynamicBody, DynamicStructNodeAd>(uL, pN, dMass_a, Xgc_a, J, fOut),
       pNode(pN)
{
}

DynamicBodyAd::~DynamicBodyAd()
{
}

void
DynamicBodyAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 12;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
DynamicBodyAd::AssJac(VariableSubMatrixHandler& WorkMat,
                      doublereal dCoef,
                      const VectorHandler& XCurr,
                      const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("DynamicBodyAd::AssJac");

     using namespace sp_grad;

     SpGradientAssVec<SpGradient>::AssJac(this,
                                          WorkMat.SetSparseGradient(),
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);

     return WorkMat;
}

void
DynamicBodyAd::AssJac(VectorHandler& JacY,
                      const VectorHandler& Y,
                      doublereal dCoef,
                      const VectorHandler& XCurr,
                      const VectorHandler& XPrimeCurr,
                      VariableSubMatrixHandler& WorkMat)
{
     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

SubVectorHandler&
DynamicBodyAd::AssRes(SubVectorHandler& WorkVec,
                      doublereal dCoef,
                      const VectorHandler& XCurr ,
                      const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("DynamicBodyAd::AssRes");

     using namespace sp_grad;

     SpGradientAssVec<doublereal>::AssRes(this,
                                          WorkVec,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

template <typename T>
void
DynamicBodyAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                      doublereal dCoef,
                      const sp_grad::SpGradientVectorHandler<T>& XCurr,
                      const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                      sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     Vec3 GravityAcceleration;

     const bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
                                              GravityAcceleration);

     const RigidBodyKinematics* const pRBK = pNode->pGetRBK();

     const integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();

     SpColVectorA<T, 3> X, V, W;
     SpMatrixA<T, 3, 3> R;

     if (pRBK) {
          pNode->GetXCurr(X, dCoef, func);
     }

     pNode->GetVCurr(V, dCoef, func);
     pNode->GetWCurr(W, dCoef, func);
     pNode->GetRCurr(R, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     if (pRBK) {
          oDofMap.GetDofStat(X);
     }

     oDofMap.GetDofStat(V);
     oDofMap.GetDofStat(W);
     oDofMap.GetDofStat(R);
     oDofMap.Reset();

     if (pRBK) {
          oDofMap.InsertDof(X);
     }

     oDofMap.InsertDof(V);
     oDofMap.InsertDof(W);
     oDofMap.InsertDof(R);
     oDofMap.InsertDone();

     const SpColVector<T, 3> STmp_sp(R * S0, oDofMap);
     const SpMatrix<T, 3, 3> JTmp_sp(SpMatrix<T, 3, 3>(R * J0, oDofMap) * Transpose(R), oDofMap);
     const SpColVector<T, 3> f1(V * -dMass - Cross(W, STmp_sp, oDofMap), oDofMap);
     const SpColVector<T, 3> f2(-Cross(STmp_sp, V, oDofMap) - SpColVector<T, 3>(JTmp_sp * W, oDofMap), oDofMap);

     WorkVec.AddItem(iFirstPositionIndex + 1, f1);
     WorkVec.AddItem(iFirstPositionIndex + 4, f2);

     if (g) {
          const SpColVector<T, 3> f3 = GravityAcceleration * dMass;
          const SpColVector<T, 3> f4(Cross(STmp_sp, GravityAcceleration), oDofMap);

          WorkVec.AddItem(iFirstPositionIndex + 7, f3);
          WorkVec.AddItem(iFirstPositionIndex + 10, f4);
     }

     if (pRBK) {
          AssVecRBK_int(pRBK, X, V, W, STmp_sp, JTmp_sp, WorkVec, oDofMap);
     }

     UpdateInertia(STmp_sp, JTmp_sp);
}

void DynamicBodyAd::UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& S,
                                  const sp_grad::SpMatrix<doublereal, 3, 3>& J) const
{
     BodyAd::UpdateInertia(S, J);

     pNode->AddInertia(dMass, STmp, JTmp);
}

StaticBodyAd::StaticBodyAd(unsigned int uL,
                           const StaticStructNodeAd* pN,
                           doublereal dMass_a,
                           const Vec3& Xgc_a,
                           const Mat3x3& J,
                           flag fOut)
     : BodyAd<StaticBody, StaticStructNodeAd>(uL, pN, dMass_a, Xgc_a, J, fOut),
       pNode(pN)
{
}


StaticBodyAd::~StaticBodyAd()
{
}

void
StaticBodyAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 12;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
StaticBodyAd::AssJac(VariableSubMatrixHandler& WorkMat,
                     doublereal dCoef,
                     const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("StaticBodyAd::AssJac");

     using namespace sp_grad;

     SpGradientAssVec<SpGradient>::AssJac(this,
                                          WorkMat.SetSparseGradient(),
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

void
StaticBodyAd::AssJac(VectorHandler& JacY,
                     const VectorHandler& Y,
                     doublereal dCoef,
                     const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr,
                     VariableSubMatrixHandler& WorkMat)
{
     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

SubVectorHandler&
StaticBodyAd::AssRes(SubVectorHandler& WorkVec,
                     doublereal dCoef,
                     const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("StaticBodyAd::AssRes");

     using namespace sp_grad;

     SpGradientAssVec<doublereal>::AssRes(this,
                                          WorkVec,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

template <typename T>
void
StaticBodyAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                     doublereal dCoef,
                     const sp_grad::SpGradientVectorHandler<T>& XCurr,
                     const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                     sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     Vec3 Acceleration(Zero3);
     const bool g = GravityOwner::bGetGravity(pNode->GetXCurr(), Acceleration);

     /* W is uninitialized because its use is conditioned by w */
     const RigidBodyKinematics* const pRBK = pNode->pGetRBK();

     SpColVectorA<T, 3> X, V, W;
     SpMatrixA<T, 3, 3> R;

     if (pRBK) {
          pNode->GetXCurr(X, dCoef, func);
          pNode->GetVCurr(V, dCoef, func);
          pNode->GetWCurr(W, dCoef, func);
     }

     pNode->GetRCurr(R, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     if (pRBK) {
          oDofMap.GetDofStat(X);
          oDofMap.GetDofStat(V);
          oDofMap.GetDofStat(W);
     }

     oDofMap.GetDofStat(R);
     oDofMap.Reset();

     if (pRBK) {
          oDofMap.InsertDof(X);
          oDofMap.InsertDof(V);
          oDofMap.InsertDof(W);
     }

     oDofMap.InsertDof(R);
     oDofMap.InsertDone();

     const SpColVector<T, 3> STmp_sp(R * S0, oDofMap);
     const SpMatrix<T, 3, 3> JTmp_sp(SpMatrix<T, 3, 3>(R * J0, oDofMap) * Transpose(R), oDofMap);

     if (g) {
          integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
          SpColVector<T, 3> FTmp = Acceleration * dMass;
          SpColVector<T, 3> MTmp(Cross(STmp_sp, Acceleration), oDofMap);

          WorkVec.AddItem(iFirstMomentumIndex + 1, FTmp);
          WorkVec.AddItem(iFirstMomentumIndex + 4, MTmp);
     }

     if (pRBK) {
          AssVecRBK_int(pRBK, X, V, W, STmp_sp, JTmp_sp, WorkVec, oDofMap);
     }

     UpdateInertia(STmp_sp, JTmp_sp);
}

ModalBodyAd::ModalBodyAd(unsigned int uL,
                         const ModalNodeAd* pNode_a,
                         doublereal dMass_a,
                         const Vec3& Xgc_a,
                         const Mat3x3& J,
                         flag fOut)
     : BodyAd<ModalBody, ModalNodeAd>(uL, pNode_a, dMass_a, Xgc_a, J, fOut),
       pNode(pNode_a)
{
}

ModalBodyAd::~ModalBodyAd(void)
{
}

void
ModalBodyAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 6;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
ModalBodyAd::AssJac(VariableSubMatrixHandler& WorkMat,
                    doublereal dCoef,
                    const VectorHandler& XCurr,
                    const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("ModalBodyAd::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);

     return WorkMat;
}

void
ModalBodyAd::AssJac(VectorHandler& JacY,
                    const VectorHandler& Y,
                    doublereal dCoef,
                    const VectorHandler& XCurr,
                    const VectorHandler& XPrimeCurr,
                    VariableSubMatrixHandler& WorkMat)
{
     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

SubVectorHandler&
ModalBodyAd::AssRes(SubVectorHandler& WorkVec,
                    doublereal dCoef,
                    const VectorHandler& XCurr,
                    const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("ModalBodyAd::AssRes");

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);
     return WorkVec;
}

template <typename T>
void
ModalBodyAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const sp_grad::SpGradientVectorHandler<T>& XCurr,
                    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                    sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     Vec3 GravityAcceleration;
     const bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
                                              GravityAcceleration);

     if (pNode->pGetRBK()) {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     SpColVector<T, 3> W(3, 3);
     SpMatrix<T, 3, 3> R(3, 3, 3);
     SpColVector<T, 3> XPP_sp(3, 1), WP_sp(3, 1);

     pNode->GetWCurr(W, dCoef, func);
     pNode->GetRCurr(R, dCoef, func);
     pNode->GetXPPCurr(XPP_sp, dCoef, func);
     pNode->GetWPCurr(WP_sp, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(W);
     oDofMap.GetDofStat(R);
     oDofMap.GetDofStat(XPP_sp);
     oDofMap.GetDofStat(WP_sp);
     oDofMap.Reset();
     oDofMap.InsertDof(W);
     oDofMap.InsertDof(R);
     oDofMap.InsertDof(XPP_sp);
     oDofMap.InsertDof(WP_sp);
     oDofMap.InsertDone();

     const SpColVector<T, 3> STmp_sp(R * S0, oDofMap);
     const SpMatrix<T, 3, 3> JTmp_sp(SpMatrix<T, 3, 3>(R * J0, oDofMap) * Transpose(R), oDofMap);

     SpColVector<T, 3> F(XPP_sp * -dMass - Cross(WP_sp, STmp_sp) - Cross(W, Cross(W, STmp_sp), oDofMap), oDofMap);
     SpColVector<T, 3> M(-Cross(STmp_sp, XPP_sp) - Cross(W, JTmp_sp * W) - SpColVector<T, 3>(JTmp_sp * WP_sp, oDofMap), oDofMap);

     if (g) {
          F += GravityAcceleration * dMass;
          M.Add(Cross(STmp_sp, GravityAcceleration), oDofMap);
     }

     const integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();

     WorkVec.AddItem(iFirstPositionIndex + 7, F);
     WorkVec.AddItem(iFirstPositionIndex + 10, M);

     UpdateInertia(STmp_sp, JTmp_sp);
}

void ModalBodyAd::UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& S,
                                const sp_grad::SpMatrix<doublereal, 3, 3>& J) const
{
     BodyAd::UpdateInertia(S, J);

     pNode->AddInertia(dMass, STmp, JTmp);
}
