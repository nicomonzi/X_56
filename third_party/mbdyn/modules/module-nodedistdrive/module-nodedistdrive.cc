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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
	Copyright (C) 2015(-2023) all rights reserved.

	The copyright of this code is transferred
	to Pierangelo Masarati and Paolo Mantegazza
	for use in the software MBDyn as described
	in the GNU Public License version 2.1
*/

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"
#endif

#include <drive.h>

#include <myassert.h>
#include <except.h>

#include <strnode.h>
#include <elem.h>
#include <mynewmem.h>
#include <dataman.h>

#include "module-nodedistdrive.h"

class NodeDistDriveCaller : public DriveCaller
{
public:
	enum OutputRefFrame {
		RFM_Node1,
		RFM_Node2
	};

	NodeDistDriveCaller(const DriveHandler* pDH,
			    const StructNode* pNode1,
			    const Vec3& o1,
			    const StructNode* pNode2,
			    const Vec3& o2,
			    const Vec3& e1,
			    OutputRefFrame eOutputRefFrame);
	virtual ~NodeDistDriveCaller(void);
	virtual DriveCaller* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	inline doublereal dGet(const doublereal& dVar) const;
	inline doublereal dGet(void) const;
	virtual bool bIsDifferentiable(void) const;
	inline virtual doublereal dGetP(void) const;
	inline virtual doublereal dGetP(const doublereal& dVar) const;
private:
	const StructNode* const pNode1;
	const Vec3 o1;
	const StructNode* const pNode2;
	const Vec3 o2;
	const Vec3 e1;
	const OutputRefFrame eOutputRefFrame;
};

class NodeRotDriveCaller : public DriveCaller
{
public:
	NodeRotDriveCaller(const DriveHandler* pDH,
			   const StructNode* pNode1,
			   const Mat3x3& e1,
			   const StructNode* pNode2,
			   const Mat3x3& e2);
	virtual ~NodeRotDriveCaller(void);
	virtual DriveCaller* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	inline doublereal dGet(const doublereal& dVar) const;
	inline doublereal dGet(void) const;
	virtual bool bIsDifferentiable(void) const;
	inline virtual doublereal dGetP(void) const;
	inline virtual doublereal dGetP(const doublereal& dVar) const;

private:
	const StructNode* const pNode1;
	const Mat3x3 e1;
	const StructNode* const pNode2;
	const Mat3x3 e2;
};


struct NodeDistDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
NodeDistDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "node distance");

	const DriveHandler* pDrvHdl = 0;

	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	if (pDM == 0) {
		silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
			<< "no DOF dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if ( !HP.IsKeyWord("node1") )
	{
		silent_cerr("node distance drive caller: keyword \"node1\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// Allow us to use dymmy nodes
	StructNode* const pNode1 = dynamic_cast<StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));

	if (!pNode1) {
		silent_cerr("A structural node is required at line " << HP.GetLineData() << "\n");
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const Vec3 o1 = HP.IsKeyWord("offset") ? HP.GetPosRel(ReferenceFrame(pNode1)) : Zero3;

	if ( !HP.IsKeyWord("node2") )
	{
		silent_cerr("node distance drive caller: keyword \"node2\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// Allow us to use dummy nodes
	StructNode* const pNode2 = dynamic_cast<StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));

	if (!pNode2) {
		silent_cerr("A structural node is required at line " << HP.GetLineData() << "\n");
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const Vec3 o2 = HP.IsKeyWord("offset") ? HP.GetPosRel(ReferenceFrame(pNode2)) : Zero3;

	NodeDistDriveCaller::OutputRefFrame eOutputRefFrame(NodeDistDriveCaller::RFM_Node2);

	if (HP.IsKeyWord("output" "reference" "frame"))
	{
		if (HP.IsKeyWord("node1")) {
			eOutputRefFrame = NodeDistDriveCaller::RFM_Node1;
		} else if (HP.IsKeyWord("node2")) {
			eOutputRefFrame = NodeDistDriveCaller::RFM_Node2;
		} else {
			silent_cerr("node distance drive caller: keyword \"node1\" or \"node2\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if ( !HP.IsKeyWord("direction") )
	{
		silent_cerr("node distance drive caller: keyword \"direction\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 e1 = HP.GetVecRel(ReferenceFrame(eOutputRefFrame == NodeDistDriveCaller::RFM_Node1 ? pNode1 : pNode2));

	e1 /= e1.Norm();

	SAFENEWWITHCONSTRUCTOR(pDC,
			       NodeDistDriveCaller,
			       NodeDistDriveCaller(pDrvHdl, pNode1, o1, pNode2, o2, e1, eOutputRefFrame));

	pDM->GetLogFile()
		<< "nodedistdrive: " << pDC->GetLabel()
		<< " (" << pDC->GetName() << ") "
		<< pNode1->GetLabel() << " "
		<< o1 << " "
		<< pNode2->GetLabel() << " "
		<< o2 << " "
		<< e1 << std::endl;

	return pDC;
}

doublereal
NodeDistDriveCaller::dGet(const doublereal& dVar) const
{
	return dGet();
}

doublereal
NodeDistDriveCaller::dGet(void) const
{
	const Vec3& X1 = pNode1->GetXCurr();
	const Mat3x3& R1 = pNode1->GetRCurr();
	const Vec3& X2 = pNode2->GetXCurr();
	const Mat3x3& R2 = pNode2->GetRCurr();
	const Mat3x3& Rout = RFM_Node1 == eOutputRefFrame ? R1 : R2;

	const doublereal dX = e1.Dot(Rout.MulTV(X2 + R2 * o2 - X1 - R1 * o1));

	return dX;
}

bool NodeDistDriveCaller::bIsDifferentiable(void) const
{
	return true;
}

doublereal NodeDistDriveCaller::dGetP(void) const
{
	const Vec3& X1 = pNode1->GetXCurr();
	const Vec3& X1Dot = pNode1->GetVCurr();
	const Mat3x3& R1 = pNode1->GetRCurr();
	const Vec3& omega1 = pNode1->GetWCurr();

	const Vec3& X2 = pNode2->GetXCurr();
	const Vec3& X2Dot = pNode2->GetVCurr();
	const Mat3x3& R2 = pNode2->GetRCurr();
	const Vec3& omega2 = pNode2->GetWCurr();

	const Mat3x3& Rout = RFM_Node1 == eOutputRefFrame ? R1 : R2;
	const Vec3& omegaout = RFM_Node1 == eOutputRefFrame ? omega1 : omega2;

	const doublereal dXP = e1.Dot(Rout.MulTV(-omegaout.Cross(X2 + R2 * o2 - X1 - R1 * o1) + X2Dot + omega2.Cross(R2 * o2) - X1Dot - omega1.Cross(R1 * o1)));

	return dXP;
}

doublereal NodeDistDriveCaller::dGetP(const doublereal& dVar) const
{
	return dGetP();
}

NodeDistDriveCaller::NodeDistDriveCaller(
		const DriveHandler* pDH,
		const StructNode* pNode1Arg,
		const Vec3& o1Arg,
		const StructNode* pNode2Arg,
		const Vec3& o2Arg,
		const Vec3& e1Arg,
		OutputRefFrame eOutputRefFrameArg)
: DriveCaller(pDH),
  pNode1(pNode1Arg), o1(o1Arg),
  pNode2(pNode2Arg), o2(o2Arg),
  e1(e1Arg),
  eOutputRefFrame(eOutputRefFrameArg)
{
	NO_OP;
};

NodeDistDriveCaller::~NodeDistDriveCaller(void)
{
	NO_OP;
}

DriveCaller*
NodeDistDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = NULL;

	SAFENEWWITHCONSTRUCTOR(pDC,
			NodeDistDriveCaller,
			NodeDistDriveCaller(
				pDrvHdl,
				pNode1, o1,
				pNode2, o2,
				e1,
				eOutputRefFrame));

	return pDC;
}

std::ostream&
NodeDistDriveCaller::Restart(std::ostream& out) const
{
	out << "node distance,"
		<< pNode1->GetLabel() << ","
		<< o1 << ","
		<< pNode2->GetLabel() << ","
		<< o2 << ","
		<< e1;

	return out;
}

struct NodeRotDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
NodeRotDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "node rotation");

	const DriveHandler* pDrvHdl = 0;

	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	if (pDM == 0) {
		silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
			<< "no DOF dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if ( !HP.IsKeyWord("node1") )
	{
		silent_cerr("node distance drive caller: keyword \"node1\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	StructNode* const pNode1 = pDM->ReadNode<StructNode, StructDispNode, Node::STRUCTURAL>(HP);

	const Mat3x3 e1 = HP.IsKeyWord("orientation") ? HP.GetRotRel(ReferenceFrame(pNode1)) : Eye3;

	if ( !HP.IsKeyWord("node2") )
	{
		silent_cerr("node distance drive caller: keyword \"node2\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	StructNode* const pNode2 = pDM->ReadNode<StructNode, StructDispNode, Node::STRUCTURAL>(HP);

	const Mat3x3 e2 = HP.IsKeyWord("orientation") ? HP.GetRotRel(ReferenceFrame(pNode2)) : Eye3;

	SAFENEWWITHCONSTRUCTOR(pDC,
			       NodeRotDriveCaller,
			       NodeRotDriveCaller(pDrvHdl, pNode1, e1, pNode2, e2));

	pDM->GetLogFile()
	    << "noderotdrive: " << pDC->GetLabel()
	    << " (" << pDC->GetName() << ") "
	    << pNode1->GetLabel() << " "
	    << e1 << " "
	    << pNode2->GetLabel() << " "
	    << e2 << std::endl;

	return pDC;
}

doublereal
NodeRotDriveCaller::dGet(const doublereal& dVar) const
{
	return dGet();
}

doublereal
NodeRotDriveCaller::dGet(void) const
{
    const Mat3x3& R1 = pNode1->GetRCurr();
    const Mat3x3& R2 = pNode2->GetRCurr();

    const Vec3 v = e2.MulTV(R2.MulTV(R1 * e1.GetCol(1)));

    return atan2(v(2), v(1));
}

bool NodeRotDriveCaller::bIsDifferentiable(void) const
{
	return true;
}

doublereal NodeRotDriveCaller::dGetP(void) const
{
    const Vec3& omega1 = pNode1->GetWCurr();
    const Vec3& omega2 = pNode2->GetWCurr();
    const Mat3x3& R2 = pNode2->GetRCurr();

    return e2.GetCol(3).Dot(R2.MulTV(omega1 - omega2));
}

doublereal NodeRotDriveCaller::dGetP(const doublereal& dVar) const
{
	return dGetP();
}

NodeRotDriveCaller::NodeRotDriveCaller(
		const DriveHandler* pDH,
		const StructNode* pNode1_a,
		const Mat3x3& e1_a,
		const StructNode* pNode2_a,
		const Mat3x3& e2_a)
: DriveCaller(pDH),
  pNode1(pNode1_a), e1(e1_a),
  pNode2(pNode2_a), e2(e2_a)
{
	NO_OP;
};

NodeRotDriveCaller::~NodeRotDriveCaller(void)
{
	NO_OP;
}

DriveCaller*
NodeRotDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = NULL;

	SAFENEWWITHCONSTRUCTOR(pDC,
			NodeRotDriveCaller,
			NodeRotDriveCaller(
				pDrvHdl,
				pNode1, e1,
				pNode2, e2));

	return pDC;
}

std::ostream&
NodeRotDriveCaller::Restart(std::ostream& out) const
{
	out << "node rotation,"
		<< pNode1->GetLabel() << ","
		<< e1 << ","
		<< pNode2->GetLabel() << ","
		<< e2;

	return out;
}

bool nodedistdrive_set()
{
	DriveCallerRead	*rf = new NodeDistDCR;

	if (!SetDriveCallerData("node" "distance", rf)) {
		delete rf;
		return false;
	}

	rf = new NodeRotDCR;

	if (!SetDriveCallerData("node" "rotation", rf)) {
		delete rf;
		return false;
	}

    return true;
}

#ifndef STATIC_MODULES

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!nodedistdrive_set()) {
		silent_cerr("nodedistdrive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

#endif // ! STATIC_MODULES
