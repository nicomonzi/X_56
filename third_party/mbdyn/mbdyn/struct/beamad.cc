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
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "beamad.h"
#include "shapefnc.h"

template <typename BeamType>
BeamAdBase<BeamType>::BeamAdBase(unsigned int uL,
                                 const StructNodeAd* pN1,
                                 const StructNodeAd* pN2,
                                 const StructNodeAd* pN3,
                                 const Vec3& F1,
                                 const Vec3& F2,
                                 const Vec3& F3,
                                 const Mat3x3& R1,
                                 const Mat3x3& R2,
                                 const Mat3x3& R3,
                                 const Mat3x3& r_I,
                                 const Mat3x3& rII,
                                 ConstitutiveLaw6D* pD_I,
                                 ConstitutiveLaw6D* pDII,
                                 OrientationDescription ood,
                                 flag fOut)
:
     BeamType(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, ood, fOut),
     pNode{pN1, pN2, pN3}
{

}

template <typename BeamType>
BeamAdBase<BeamType>::~BeamAdBase()
{
}

template <typename BeamType>
void BeamAdBase<BeamType>::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 18;
     *piNumCols = 0;
}

template <typename BeamType>
void
BeamAdBase<BeamType>::AddInternalForcesAD(sp_grad::SpColVector<doublereal, 6>& AzLoc_a, unsigned int iSez)
{
     // Not in use so far, but may be overridden by derived classes.
}

template <typename BeamType>
void
BeamAdBase<BeamType>::AddInternalForcesAD(sp_grad::SpColVector<sp_grad::SpGradient, 6>& AzLoc_a, unsigned int iSez)
{
}

template <typename BeamType>
void
BeamAdBase<BeamType>::AddInternalForcesAD(sp_grad::SpColVector<sp_grad::GpGradProd, 6>& AzLoc_a, unsigned int iSez)
{
}

template <typename BeamType>
template <typename T>
void
BeamAdBase<BeamType>::InterpState(const sp_grad::SpColVector<T, 3>& v1,
                                  const sp_grad::SpColVector<T, 3>& v2,
                                  const sp_grad::SpColVector<T, 3>& v3,
                                  sp_grad::SpColVector<T, 3>& p,
                                  Section Sec,
                                  const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= 3; ++i) {
          oDofMap.MapAssign(p(i), v1(i) * dN3[Sec][0] + v2(i) * dN3[Sec][1] + v3(i) * dN3[Sec][2]);
     }
}

template <typename BeamType>
template <typename T>
void
BeamAdBase<BeamType>::InterpDeriv(const sp_grad::SpColVector<T, 3>& v1,
                                  const sp_grad::SpColVector<T, 3>& v2,
                                  const sp_grad::SpColVector<T, 3>& v3,
                                  sp_grad::SpColVector<T, 3>& g_a,
                                  Section Sec,
                                  const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= 3; ++i) {
          oDofMap.MapAssign(g_a(i), (v1(i) * dN3P[Sec][0] + v2(i) * dN3P[Sec][1] + v3(i) * dN3P[Sec][2]) * dsdxi[Sec]);
     }
}

template <typename BeamType>
template <typename T>
void
BeamAdBase<BeamType>::AssReactionForce(sp_grad::SpGradientAssVec<T>& WorkVec,
                                       const std::array<sp_grad::SpColVectorA<T, 3>, NUMSEZ>& p_a,
                                       const std::array<sp_grad::SpColVectorA<T, 6>, NUMSEZ>& Az_a,
                                       const std::array<sp_grad::SpColVectorA<T, 3>, NUMNODES>& X,
                                       const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) const
{
     using namespace sp_grad;

     const index_type iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
     const index_type iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
     const index_type iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();

     const SpColVector<T, 3> F_I = SubColVector<1, 1, 3>(Az_a[S_I]);

     DEBUGCERR("BeamAd(" << GetLabel() << "): F_I=" << F_I << "\n");

     WorkVec.AddItem(iNode1FirstMomIndex + 1, F_I);

     const SpColVector<T, 3> M_I(Cross(p_a[S_I] - X[NODE1], SubColVector<1, 1, 3>(Az_a[S_I]), oDofMap) + SubColVector<4, 1, 3>(Az_a[S_I]), oDofMap);

     DEBUGCERR("BeamAd(" << GetLabel() << "): M_I=" << M_I << "\n");

     WorkVec.AddItem(iNode1FirstMomIndex + 4, M_I);

     const SpColVector<T, 3> F_II(SubColVector<1, 1, 3>(Az_a[SII]) - SubColVector<1, 1, 3>(Az_a[S_I]), oDofMap);
     const SpColVector<T, 3> M_II(SubColVector<4, 1, 3>(Az_a[SII]) - SubColVector<4, 1, 3>(Az_a[S_I])
                                  + Cross(p_a[SII] - X[NODE2], SubColVector<1, 1, 3>(Az_a[SII]), oDofMap)
                                  - Cross(p_a[S_I] - X[NODE2], SubColVector<1, 1, 3>(Az_a[S_I]), oDofMap), oDofMap);

     DEBUGCERR("BeamAd(" << GetLabel() << "): F_II=" << F_II << "\n");

     WorkVec.AddItem(iNode2FirstMomIndex + 1, F_II);

     DEBUGCERR("BeamAd(" << GetLabel() << "): M_II=" << M_II << "\n");

     WorkVec.AddItem(iNode2FirstMomIndex + 4, M_II);

     const SpColVector<T, 3> F_III = -SubColVector<1, 1, 3>(Az_a[SII]);
     const SpColVector<T, 3> M_III(Cross(SubColVector<1, 1, 3>(Az_a[SII]), p_a[SII] - X[NODE3], oDofMap) - SubColVector<4, 1, 3>(Az_a[SII]), oDofMap);

     WorkVec.AddItem(iNode3FirstMomIndex + 1, F_III);
     WorkVec.AddItem(iNode3FirstMomIndex + 4, M_III);
}

BeamAd::BeamAd(unsigned int uL,
               const StructNodeAd* pN1,
               const StructNodeAd* pN2,
               const StructNodeAd* pN3,
               const Vec3& F1,
               const Vec3& F2,
               const Vec3& F3,
               const Mat3x3& R1,
               const Mat3x3& R2,
               const Mat3x3& R3,
               const Mat3x3& r_I,
               const Mat3x3& rII,
               ConstitutiveLaw6D* pD_I,
               ConstitutiveLaw6D* pDII,
               OrientationDescription ood,
               flag fOut)
:
     BeamAdBase<Beam>(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, ood, fOut)
// pNode{pN1, pN2, pN3}
{

}

BeamAd::~BeamAd()
{
}

void BeamAd::UpdateState(const std::array<sp_grad::SpMatrixA<doublereal, 3, 3>, NUMSEZ>& RTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& pTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& gTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& LTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefLocTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzLocTmp)
{
     using namespace sp_grad;

     for (index_type i = 0; i < NUMSEZ; ++i) {
          R[i] = RTmp[i];
          p[i] = pTmp[i];
          g[i] = gTmp[i];
          L[i] = LTmp[i];
          DefLoc[i] = DefLocTmp[i];
          Az[i] = AzTmp[i];
          AzLoc[i] = AzLocTmp[i];
     }
}

template <typename T>
inline void
BeamAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
               doublereal dCoef,
               const sp_grad::SpGradientVectorHandler<T>& XCurr,
               const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
               enum sp_grad::SpFunctionCall func)
{
     UnivAssRes(WorkVec, dCoef, XCurr, func);
}

SubVectorHandler&
BeamAd::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
     sp_grad::SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                          WorkVec,
                                                          XCurr,
                                                          sp_grad::SpFunctionCall::INITIAL_ASS_RES);

     return WorkVec;
}

VariableSubMatrixHandler&
BeamAd::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr)
{
     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::InitialAssJac(this,
                                                                   WorkMat.SetSparseGradient(),
                                                                   XCurr,
                                                                   sp_grad::SpFunctionCall::INITIAL_ASS_JAC);

     return WorkMat;
}

template <typename T>
inline void
BeamAd::InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                      const sp_grad::SpGradientVectorHandler<T>& XCurr,
                      sp_grad::SpFunctionCall func)
{
     UnivAssRes(WorkVec, 1., XCurr, func);
}

template <typename T>
inline void
BeamAd::UnivAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                   doublereal dCoef,
                   const sp_grad::SpGradientVectorHandler<T>& XCurr,
                   enum sp_grad::SpFunctionCall func)
{
     DEBUGCOUTFNAME("BeamAd::UnivAssRes");

     using namespace sp_grad;

     DEBUGCOUT("dCoef=" << dCoef << "\n");

     std::array<SpColVectorA<T, 3>, NUMNODES> X, gNod;

     SpGradExpDofMapHelper<T> oDofMap;

     for (unsigned int i = 0; i < NUMNODES; i++) {
          pNode[i]->GetgCurr(gNod[i], dCoef, func);
          pNode[i]->GetXCurr(X[i], dCoef, func);

          oDofMap.GetDofStat(X[i]);
          oDofMap.GetDofStat(gNod[i]);
     }

     oDofMap.Reset();

     for (unsigned int i = 0; i < NUMNODES; ++i) {
          oDofMap.InsertDof(X[i]);
          oDofMap.InsertDof(gNod[i]);
     }

     oDofMap.InsertDone();

     SpMatrixA<T, 3, 3> RNod;
     std::array<SpColVectorA<T, 3>, NUMNODES> xTmp;

     for (unsigned int i = 0; i < NUMNODES; i++) {
          pNode[i]->GetRCurr(RNod, dCoef, func); // No need to insert into oDofMap because it will depend only on gNod
          xTmp[i].MapAssign(X[i] + RNod * f[i], oDofMap);

          DEBUGCOUT("f[" << i << "]=" << f[i] << "\n");
          DEBUGCOUT("gNod[" << i << "]=" << gNod[i] << "\n");
          DEBUGCOUT("RNod[" << i << "]=" << RNod << "\n");
          DEBUGCOUT("X[" << i << "]=" << X[i] << "\n");
          DEBUGCOUT("xTmp[" << i << "]=" << xTmp[i] << "\n");
     }

     std::array<SpMatrixA<T, 3, 3>, NUMSEZ> RDelta, R_sp;
     std::array<SpColVectorA<T, 3>, NUMSEZ> gGrad, p_sp, g_sp, L_sp;
     std::array<SpColVectorA<T, 6>, NUMSEZ> DefLoc_sp, Az_sp, AzLoc_sp;

     DEBUGCOUT("beam3(" << GetLabel() << ")\n");
     DEBUGCOUT("Beam::AssRes bFirstRes = " << bFirstRes << std::endl);

     /* Aggiorna le grandezze della trave nei punti di valutazione */
     for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {

          /* Posizione */
          InterpState(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], p_sp[iSez], Beam::Section(iSez), oDofMap);

          /* Matrici di rotazione */
          InterpState(gNod[NODE1], gNod[NODE2], gNod[NODE3], g_sp[iSez], Beam::Section(iSez), oDofMap);
          MatRVec(g_sp[iSez], RDelta[iSez], oDofMap);
          R_sp[iSez].MapAssign(RDelta[iSez] * RRef[iSez], oDofMap);

          /* Derivate della posizione */
          InterpDeriv(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], L_sp[iSez], Beam::Section(iSez), oDofMap);

          /* Derivate dei parametri di rotazione */
          InterpDeriv(gNod[NODE1], gNod[NODE2], gNod[NODE3], gGrad[iSez], Beam::Section(iSez), oDofMap);

          /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
          const SpColVector<T, 3> GgGrad(MatGVec(g_sp[iSez], oDofMap) * gGrad[iSez], oDofMap);

          for (index_type i = 1; i <= 3; ++i) {
               DefLoc_sp[iSez](i) = Dot(R_sp[iSez].GetCol(i), L_sp[iSez], oDofMap) - L0[iSez](i);
               DefLoc_sp[iSez](i + 3) = Dot(R_sp[iSez].GetCol(i), GgGrad, oDofMap) + DefLocRef[iSez](i + 3);
          }

          DEBUGCERR("BeamAd(" << GetLabel() << "): DefLoc[" << iSez << "]=" << DefLoc_sp[iSez] << "\n");

          /* Calcola le azioni interne */
          pD[iSez]->Update(DefLoc_sp[iSez], AzLoc_sp[iSez], oDofMap);

          /* corregge le azioni interne locali (piezo, ecc) */
          AddInternalForcesAD(AzLoc_sp[iSez], iSez);

          /* Porta le azioni interne nel sistema globale */
          for (integer i = 1; i <= 3; ++i) {
               Az_sp[iSez](i) = Dot(Transpose(R_sp[iSez].GetRow(i)), SubColVector<1, 1, 3>(AzLoc_sp[iSez]), oDofMap);
               Az_sp[iSez](i + 3) = Dot(Transpose(R_sp[iSez].GetRow(i)), SubColVector<4, 1, 3>(AzLoc_sp[iSez]), oDofMap);
          }

          DEBUGCOUT("p[" << iSez << "]=" << p_sp[iSez] << std::endl);
          DEBUGCOUT("g[" << iSez << "]=" << g_sp[iSez] << std::endl);
          DEBUGCOUT("RDelta[" << iSez << "]=" << RDelta[iSez] << std::endl);
          DEBUGCOUT("RPrev[" << iSez << "]=" << RPrev[iSez] << std::endl);
          DEBUGCOUT("RRef[" << iSez << "]=" << RRef[iSez] << std::endl);
          DEBUGCOUT("R[" << iSez << "]=" << R_sp[iSez] << std::endl);
          DEBUGCOUT("L[" << iSez << "]=" << L_sp[iSez] << std::endl);
          DEBUGCOUT("L0[" << iSez << "]=" << L0[iSez] << std::endl);
          DEBUGCOUT("DefLoc[" << iSez << "]=" << DefLoc_sp[iSez] << std::endl);
          DEBUGCOUT("DefLocRef[" << iSez << "]=" << DefLocRef[iSez] << std::endl);
          DEBUGCOUT("Az[" << iSez << "]=" << Az_sp[iSez] << std::endl);
          DEBUGCOUT("AzLoc[" << iSez << "]=" << AzLoc_sp[iSez] << std::endl);
     }

     AssReactionForce(WorkVec, p_sp, Az_sp, X, oDofMap);

     UpdateState(R_sp, p_sp, g_sp, L_sp, DefLoc_sp, Az_sp, AzLoc_sp);

     bFirstRes = false;
}

VariableSubMatrixHandler&
BeamAd::AssJac(VariableSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("BeamAd::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

void
BeamAd::AssJac(VectorHandler& JacY,
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
BeamAd::AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("BeamAd::AssRes");

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

ViscoElasticBeamAd::ViscoElasticBeamAd(unsigned int uL,
                                       const StructNodeAd* pN1,
                                       const StructNodeAd* pN2,
                                       const StructNodeAd* pN3,
                                       const Vec3& F1,
                                       const Vec3& F2,
                                       const Vec3& F3,
                                       const Mat3x3& R1,
                                       const Mat3x3& R2,
                                       const Mat3x3& R3,
                                       const Mat3x3& r_I,
                                       const Mat3x3& rII,
                                       ConstitutiveLaw6D* pD_I,
                                       ConstitutiveLaw6D* pDII,
                                       OrientationDescription ood,
                                       flag fOut)
:
     BeamAdBase<ViscoElasticBeam>(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, ood, fOut)
{
}

ViscoElasticBeamAd::~ViscoElasticBeamAd()
{
}

inline void
ViscoElasticBeamAd::UpdateState(const std::array<sp_grad::SpMatrixA<doublereal, 3, 3>, NUMSEZ>& RTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& pTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& gTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& gPrimeTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& OmegaTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& LTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& LPrimeTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefLocTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefPrimeLocTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzLocTmp)
{
     using namespace sp_grad;

     for (index_type i = 0; i < NUMSEZ; ++i) {
          R[i] = RTmp[i];
          p[i] = pTmp[i];
          g[i] = gTmp[i];
          gPrime[i] = gPrimeTmp[i];
          L[i] = LTmp[i];
          LPrime[i] = LPrimeTmp[i];
          DefLoc[i] = DefLocTmp[i];
          DefPrimeLoc[i] = DefPrimeLocTmp[i];
          Az[i] = AzTmp[i];
          AzLoc[i] = AzLocTmp[i];
     }
}

template <typename T>
inline void
ViscoElasticBeamAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                           doublereal dCoef,
                           const sp_grad::SpGradientVectorHandler<T>& XCurr,
                           const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                           enum sp_grad::SpFunctionCall func)
{
     UnivAssRes(WorkVec, dCoef, XCurr, func);
}

SubVectorHandler&
ViscoElasticBeamAd::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
     sp_grad::SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                          WorkVec,
                                                          XCurr,
                                                          sp_grad::SpFunctionCall::INITIAL_ASS_RES);

     return WorkVec;
}

VariableSubMatrixHandler&
ViscoElasticBeamAd::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                  const VectorHandler& XCurr)
{
     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::InitialAssJac(this,
                                                                   WorkMat.SetSparseGradient(),
                                                                   XCurr,
                                                                   sp_grad::SpFunctionCall::INITIAL_ASS_JAC);

     return WorkMat;
}

template <typename T>
inline void
ViscoElasticBeamAd::InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                  const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                  sp_grad::SpFunctionCall func)
{
     UnivAssRes(WorkVec, 1., XCurr, func);
}

template <typename T>
inline void
ViscoElasticBeamAd::UnivAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                               doublereal dCoef,
                               const sp_grad::SpGradientVectorHandler<T>& XCurr,
                               enum sp_grad::SpFunctionCall func)
{
     DEBUGCOUTFNAME("ViscoElasticBeamAd::UnivAssRes");
     /* Riceve il vettore gia' dimensionato e con gli indici a posto
      * per scrivere il residuo delle equazioni di equilibrio dei tre nodi */

     /* Per la trattazione teorica, il riferimento e' il file ul-travi.tex
      * (ora e' superato) */

     using namespace sp_grad;

     std::array<SpColVectorA<T, 3>, NUMNODES> gNod, gPrimeNod;
     std::array<SpColVectorA<T, 3>, NUMNODES> XNod, XPrimeNod;

     SpGradExpDofMapHelper<T> oDofMap;

     for (unsigned int i = 0; i < NUMNODES; i++) {
          pNode[i]->GetgCurr(gNod[i], dCoef, func);
          pNode[i]->GetgPCurr(gPrimeNod[i], dCoef, func);
          pNode[i]->GetXCurr(XNod[i], dCoef, func);
          pNode[i]->GetVCurr(XPrimeNod[i], dCoef, func);

          oDofMap.GetDofStat(XNod[i]);
          oDofMap.GetDofStat(gNod[i]);
          oDofMap.GetDofStat(XPrimeNod[i]);
          oDofMap.GetDofStat(gPrimeNod[i]);
     }

     oDofMap.Reset();

     for (unsigned i = 0; i < NUMNODES; ++i) {
          oDofMap.InsertDof(XNod[i]);
          oDofMap.InsertDof(gNod[i]);
          oDofMap.InsertDof(XPrimeNod[i]);
          oDofMap.InsertDof(gPrimeNod[i]);
     }

     oDofMap.InsertDone();

     std::array<SpColVectorA<T, 3>, NUMNODES> xTmp, xPrimeTmp;
     SpMatrixA<T, 3, 3> RNod;
     SpColVectorA<T, 3> WNod;

     for (unsigned int i = 0; i < NUMNODES; i++) {
          pNode[i]->GetRCurr(RNod, dCoef, func); // No need to insert into oDofMap because it will depend only on gNod
          pNode[i]->GetWCurr(WNod, dCoef, func);

          const SpColVector<T, 3> fTmp(RNod * f[i], oDofMap);

          xTmp[i].MapAssign(XNod[i] + fTmp, oDofMap);
          xPrimeTmp[i].MapAssign(XPrimeNod[i] + Cross(WNod, fTmp), oDofMap);
     }

     std::array<SpMatrixA<T, 3, 3>, NUMSEZ> R_sp, RDelta;
     std::array<SpColVectorA<T, 3>, NUMSEZ> p_sp, g_sp, gGrad, gPrime_sp, gPrimeGrad, Omega_sp, L_sp, LPrime_sp;
     std::array<SpColVectorA<T, 6>, NUMSEZ> DefLoc_sp, DefPrimeLoc_sp, Az_sp, AzLoc_sp;

     /* Aggiorna le grandezze della trave nei punti di valutazione */
     for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {

          /* Posizione */
          InterpState(xTmp[NODE1],
                      xTmp[NODE2],
                      xTmp[NODE3],
                      p_sp[iSez],
                      Beam::Section(iSez),
                      oDofMap);

          /* Matrici di rotazione */
          InterpState(gNod[NODE1],
                      gNod[NODE2],
                      gNod[NODE3],
                      g_sp[iSez],
                      Beam::Section(iSez),
                      oDofMap);

          const SpMatrix<T, 3, 3> G = MatGVec(g_sp[iSez], oDofMap);

          MatRVec(g_sp[iSez], RDelta[iSez], oDofMap);
          R_sp[iSez].MapAssign(RDelta[iSez] * RRef[iSez], oDofMap);

          /* Velocita' angolare della sezione */
          InterpState(gPrimeNod[NODE1],
                      gPrimeNod[NODE2],
                      gPrimeNod[NODE3],
                      gPrime_sp[iSez],
                      Beam::Section(iSez),
                      oDofMap);

          Omega_sp[iSez].MapAssign(G * gPrime_sp[iSez]
                                + RDelta[iSez] * OmegaRef[iSez], oDofMap);

          /* rate of MatG */
          const T dtmp0 = Dot(g_sp[iSez], g_sp[iSez], oDofMap);
          const T dtmp1 = 4. + dtmp0;
          const T dtmp2 = -4. / (dtmp1 * dtmp1);
          const T dtmp3 = 2. / dtmp1;

          const SpColVector<T, 3> GPrimeg((gPrime_sp[iSez] * dtmp0 + g_sp[iSez] * Dot(gPrime_sp[iSez], g_sp[iSez], oDofMap)) * dtmp2
                                          + Cross(gPrime_sp[iSez], g_sp[iSez]) * dtmp3, oDofMap);

          /* Derivate della posizione */
          InterpDeriv(xTmp[NODE1],
                      xTmp[NODE2],
                      xTmp[NODE3],
                      L_sp[iSez],
                      Beam::Section(iSez),
                      oDofMap);

          /* Derivate della velocita' */
          InterpDeriv(xPrimeTmp[NODE1],
                      xPrimeTmp[NODE2],
                      xPrimeTmp[NODE3],
                      LPrime_sp[iSez],
                      Beam::Section(iSez),
                      oDofMap);

          /* Derivate dei parametri di rotazione */
          InterpDeriv(gNod[NODE1],
                      gNod[NODE2],
                      gNod[NODE3],
                      gGrad[iSez],
                      Beam::Section(iSez),
                      oDofMap);

          /* Derivate delle derivate spaziali dei parametri di rotazione */
          InterpDeriv(gPrimeNod[NODE1],
                      gPrimeNod[NODE2],
                      gPrimeNod[NODE3],
                      gPrimeGrad[iSez],
                      Beam::Section(iSez),
                      oDofMap);

          /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
          const SpColVector<T, 3> GgGrad(G * gGrad[iSez], oDofMap);

          for (index_type i = 1; i <= 3; ++i) {
               DefLoc_sp[iSez](i) = Dot(R_sp[iSez].GetCol(i), L_sp[iSez], oDofMap) - L0[iSez](i);
               DefLoc_sp[iSez](i + 3) = Dot(R_sp[iSez].GetCol(i), GgGrad, oDofMap) + DefLocRef[iSez](i + 3);
          }

          DEBUGCERR("ViscoElasticBeamAd(" << GetLabel() << "): DefLoc[" << iSez << "]=" << DefLoc_sp[iSez] << "\n");

          /* Calcola le velocita' di deformazione nel sistema locale nei punti di valutazione */
          const SpColVector<T, 3> DL1(LPrime_sp[iSez] + Cross(L_sp[iSez], Omega_sp[iSez]), oDofMap);
          const SpColVector<T, 3> DL2(G * gPrimeGrad[iSez] + GPrimeg + Cross(GgGrad, Omega_sp[iSez]), oDofMap);

          for (index_type i = 1; i <= 3; ++i) {
               DefPrimeLoc_sp[iSez](i) = Dot(R_sp[iSez].GetCol(i), DL1, oDofMap);
               DefPrimeLoc_sp[iSez](i + 3) = Dot(R_sp[iSez].GetCol(i), DL2, oDofMap) + DefPrimeLocRef[iSez](i + 3);
          }

          DEBUGCERR("ViscoElasticBeamAd(" << GetLabel() << "): DefPrimeLoc[" << iSez << "]=" << DefPrimeLoc_sp[iSez] << "\n");

          /* Calcola le azioni interne */
          pD[iSez]->Update(DefLoc_sp[iSez], DefPrimeLoc_sp[iSez], AzLoc_sp[iSez], oDofMap);

          /* corregge le azioni interne locali (piezo, ecc) */
          AddInternalForcesAD(AzLoc_sp[iSez], iSez);

          /* Porta le azioni interne nel sistema globale */
          for (index_type i = 1; i <= 3; ++i) {
               Az_sp[iSez](i) = Dot(Transpose(R_sp[iSez].GetRow(i)), SubColVector<1, 1, 3>(AzLoc_sp[iSez]), oDofMap);
               Az_sp[iSez](i + 3) = Dot(Transpose(R_sp[iSez].GetRow(i)), SubColVector<4, 1, 3>(AzLoc_sp[iSez]), oDofMap);
          }
     }

     AssReactionForce(WorkVec, p_sp, Az_sp, XNod, oDofMap);

     UpdateState(R_sp, p_sp, g_sp, gPrime_sp, Omega_sp, L_sp, LPrime_sp, DefLoc_sp, DefPrimeLoc_sp, Az_sp, AzLoc_sp);

     bFirstRes = false;
}

SubVectorHandler&
ViscoElasticBeamAd::AssRes(SubVectorHandler& WorkVec,
                           doublereal dCoef,
                           const VectorHandler& XCurr,
                           const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("ViscoElasticBeamAd::AssRes");

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

VariableSubMatrixHandler&
ViscoElasticBeamAd::AssJac(VariableSubMatrixHandler& WorkMat,
                           doublereal dCoef,
                           const VectorHandler& XCurr,
                           const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("ViscoElasticBeamAd::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);

     return WorkMat;
}

void
ViscoElasticBeamAd::AssJac(VectorHandler& JacY,
                           const VectorHandler& Y,
                           doublereal dCoef,
                           const VectorHandler& XCurr,
                           const VectorHandler& XPrimeCurr,
                           VariableSubMatrixHandler& WorkMat)
{
     DEBUGCOUTFNAME("ViscoElasticBeamAd::AssJac");

     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

void
ViscoElasticBeamAd::AfterConvergence(const VectorHandler& X,
                                     const VectorHandler& XP)
{
     for (unsigned i = 0; i < NUMSEZ; i++) {
          RPrev[i] = R[i];
          DefLocPrev[i] = DefLoc[i];
          pD[i]->AfterConvergence(DefLoc[i], DefPrimeLoc[i]);
     }
}
