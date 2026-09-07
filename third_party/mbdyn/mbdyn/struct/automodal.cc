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

#include "automodal.h"

AutomaticModalElem::AutomaticModalElem(const ModalNode* pN)
     :Elem(pN->GetLabel(), pN->fToBeOutput()), pModalNode(pN)
{
}

AutomaticModalElem::~AutomaticModalElem()
{
}

Elem::Type AutomaticModalElem::GetElemType() const
{
     return Elem::AUTOMATICSTRUCTURAL;
}

std::ostream& AutomaticModalElem::Restart(std::ostream& out) const
{
     return out;
}

void AutomaticModalElem::Restart(RestartData& oData, RestartData::RestartAction eAction)
{
}

void AutomaticModalElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 6;
     *piNumCols = 12;
}

VariableSubMatrixHandler&
AutomaticModalElem::AssJac(VariableSubMatrixHandler& WorkMat,
                           doublereal dCoef,
                           const VectorHandler& XCurr,
                           const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("AutomaticModalElem::AssJac");

     FullSubMatrixHandler& WM = WorkMat.SetFull();

     integer iNumRows;
     integer iNumCols;

     WorkSpaceDim(&iNumRows, &iNumCols);
     WM.ResizeReset(iNumRows, iNumCols);

     const integer iRigidIndex = pModalNode->iGetFirstIndex();

     for (unsigned int iCnt = 1; iCnt <= 6; ++iCnt) {
          WM.PutRowIndex(iCnt, iRigidIndex + iCnt);
     }

     for (unsigned int iCnt = 1; iCnt <= 12; ++iCnt) {
          WM.PutColIndex(iCnt, iRigidIndex + iCnt);
     }

     /* completa lo Jacobiano con l'aggiunta delle equazioni {xP} = {v}
        {gP} - [wr/\]{g} = {w} */
     for (int iCnt = 1; iCnt <= 6; iCnt++) {
          WM.IncCoef(iCnt, iCnt, 1.);
          WM.DecCoef(iCnt, 6 + iCnt, dCoef);
     }

     const Vec3& wr = pModalNode->GetWRef();

     WM.Sub(3 + 1, 3 + 1, Mat3x3(MatCross, wr*dCoef));

     return WorkMat;
}

SubVectorHandler&
AutomaticModalElem::AssRes(SubVectorHandler& WorkVec,
                           doublereal dCoef,
                           const VectorHandler& XCurr,
                           const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("AutomaticModalElem::AssRes");

     integer iNumRows;
     integer iNumCols;

     WorkSpaceDim(&iNumRows, &iNumCols);

     WorkVec.ResizeReset(iNumRows);

     integer iRigidIndex = pModalNode->iGetFirstIndex();

     for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
          WorkVec.PutRowIndex(iCnt, iRigidIndex + iCnt);
     }

     const Vec3& xP = pModalNode->GetVCurr();
     const Vec3& g = pModalNode->GetgCurr();
     const Vec3& gP = pModalNode->GetgPCurr();
     const Vec3& wr = pModalNode->GetWRef();

     const Vec3 v(XCurr, iRigidIndex + 6 + 1);
     const Vec3 w(XCurr, iRigidIndex + 9 + 1);

     /* equazioni per abbassare di grado il sistema */
     WorkVec.Add(1, v - xP);
     WorkVec.Add(3 + 1, w - Mat3x3(CGR_Rot::MatG, g)*gP
                 - Mat3x3(CGR_Rot::MatR, g)*wr);

     return WorkVec;
}

void AutomaticModalElem::OutputPrepare(OutputHandler &OH)
{
}

void AutomaticModalElem::Output(OutputHandler& OH) const
{
}

void AutomaticModalElem::SetValue(DataManager *pDM,
                                  VectorHandler& X, VectorHandler& XP,
                                  SimulationEntity::Hints *ph)
{
     const integer iRigidIndex = pModalNode->iGetFirstIndex();

     X.Put(iRigidIndex + 6 + 1, pModalNode->GetVCurr());
     X.Put(iRigidIndex + 9 + 1, pModalNode->GetWCurr());
     XP.Put(iRigidIndex + 6 + 1, pModalNode->GetXPPCurr());
     XP.Put(iRigidIndex + 9 + 1, pModalNode->GetWPCurr());
}

void AutomaticModalElem::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
     connectedNodes.clear();
     connectedNodes.push_back(pModalNode);
}
