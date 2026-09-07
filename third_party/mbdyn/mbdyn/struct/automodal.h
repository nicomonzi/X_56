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

#ifndef AUTOMODAL_H
#define AUTOMODAL_H

#include "elem.h"
#include "strnode.h"
#include "matvec3.h"

class AutomaticModalElem : public Elem {
public:
     explicit AutomaticModalElem(const ModalNode* pN);

     virtual Elem::Type GetElemType() const override;

     virtual ~AutomaticModalElem();

     virtual std::ostream& Restart(std::ostream& out) const override;
     virtual void Restart(RestartData& oData, RestartData::RestartAction eAction) override;

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

     virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     virtual SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     void OutputPrepare(OutputHandler &OH) override;
     virtual void Output(OutputHandler& OH) const override;

     virtual void SetValue(DataManager *pDM,
                           VectorHandler& X, VectorHandler& XP,
                           SimulationEntity::Hints *ph = 0) override;

     virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const override;
private:
     const ModalNode* const pModalNode;
};

#endif
