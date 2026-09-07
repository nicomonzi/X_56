/*
 * Copyright (C) 2006
 * Marco Morandini
 *
 */

/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cmath>
#include <limits>

#include "loadable.h"
#include "ScalarFunctions.h"
#include "friction.h"

// #define _GNU_SOURCE 1
// #include <fenv.h>
// static void __attribute__ ((constructor))
// trapfpe ()
// {
//   /* Enable some exceptions.  At startup all exceptions are masked.  */
// 
//   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
// }


/*
 * Usage:
 *
 *	loadable: <label>, <module name>, help [ , <module data> ] ;
 */

struct module_wheellugre {
	/*
	 * Ruota
	 */
	StructNode *pWheel;

	StructNode *pRim;

	/*
	 * Direzione asse ruota
	 */
	Vec3 WheelAxle;
	
	/*
	 * Terreno
	 */
	StructNode *pGround;

	/*
	 * Posizione e orientazione del terreno
	 */
	Vec3 GroundPosition;
	Vec3 GroundDirection;

	/*
	 * Geometria ruota
	 */
	doublereal dRadius;

	doublereal dRNP;		/* R+nG'*pG */

	doublereal C_cz;

	/*
	 * Attrito
	 */
	bool fSlip;
	doublereal Preload;
	const BasicScalarFunction *re_function;
	const BasicScalarFunction *F_function;
	doublereal dHystVRef;

	doublereal kb;

	DriveCaller *pMuY0;
	DriveCaller *pMuY1;

	/*
	 * Output
	 */
	Vec3 F;
	Vec3 Fdamp;
	Vec3 M;
	doublereal dInstRadius;
	doublereal dDeltaL;
	doublereal dVn;
	doublereal dSr;
	doublereal dAlpha;
	doublereal dAlphaThreshold;
	doublereal dMuX;
	doublereal dMuY;
	doublereal dVa;
	doublereal dVc;
	doublereal V_c_sx;
	
	doublereal sigma_c;
	doublereal kappa0;
	doublereal kappa;
	doublereal L;
	std::vector<doublereal> dFn;
	std::vector<doublereal> v;
	std::vector<doublereal> params;
	doublereal F_cx;
	doublereal V_cr;
	doublereal V_cx;
	doublereal re;
	doublereal a;
	doublereal mu;
	Vec3 fwd;

	BasicFriction * fc;
	
	doublereal zeta, zetap;
	
};

/* funzioni di default */
static void*
read(LoadableElem* pEl,
		DataManager* pDM,
		MBDynParser& HP)
{
	DEBUGCOUTFNAME("read");
	
	/* allocation of user-defined struct */
	module_wheellugre* p = NULL;
	SAFENEW(p, module_wheellugre);
	p->dFn.resize(1);
	p->v.resize(2);
	p->params.resize(1);

	/*
	 * help
	 */
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	wheellugre							\n"
"Author: 	Marco Morandini <morandini@aero.polimi.it>			\n"
"									\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it				\n"
"									\n"
"	All rights reserved						\n"
"									\n"
"Connects 2 structural nodes:						\n"
"     -	Wheel								\n"
"     -	Ground								\n"
"									\n"
"Note: 									\n"
"     -	The Axle and the Wheel structural nodes must be connected 	\n"
"	by a joint that allows relative rotations only about 		\n"
"	one axis (the axle)						\n"
"     -	The center of the wheel is assumed coincident with 		\n"
"	the position of the wheel structural node			\n"
"     -	The Ground structural node supports a plane defined		\n"
"	a point and a direction orthogonal to the plane (future 	\n"
"	versions might use an arbitrary, deformable surface)		\n"
"     -	The forces are applied at the \"contact point\", that 		\n"
"	is defined according to geometrical properties 			\n"
"	of the system and according to the relative position 		\n"
"	and orientation of the Wheel and Ground structural nodes	\n"
"									\n"
"     -	Input:								\n"
"		<wheel structural node label> ,				\n"
"		<wheel axle direction> ,				\n"
"		<ground structural node label> ,			\n"
"		<reference point position of the ground plane> ,	\n"
"		<direction orthogonal to the ground plane> ,		\n"
"		<wheel radius> ,					\n"
"		<vertical reaction scalar function> ,			\n"
"		[ slip ,						\n"
"		[ preload , <vertical preload> ]			\n"
"		<longitudinal friction law> ,				\n"
"		<effective radius scalar function> ,			\n"
"		<lateral friction coefficient drive for s.r.=0>		\n"
"		<lateral friction coefficient drive for s.r.=1>		\n"
"		[ , threshold , <slip angle velocity threshold> ] ]	\n"
"		[ , damping , <wheel rim structural node label> , 	\n"
"			<tread damping stiffness> ]			\n"
"									\n"
"     -	Output:								\n"
"		1)	element label					\n"
"		2-4)	tire force in global reference frame		\n"
"		5-7)	tire couple in global reference frame		\n"
"		8)	effective radius				\n"
"		9)	tire radial deformation				\n"
"		10)	tire radial deformation velocity		\n"
"		11)	slip ratio					\n"
"		12)	slip angle					\n"
"		13)	longitudinal friction coefficient		\n"
"		14)	lateral friction coefficient			\n"
"		15)	axis relative tangential velocity		\n"
"		16)	point of contact relative tangential velocity	\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	/*
	 * leggo la ruota
	 */
	p->pWheel = pDM->ReadNode<StructNode, Node::STRUCTURAL>(HP);

	/*
	 * leggo l'orientazione dell'asse ruota nel sistema locale
	 */
	ReferenceFrame RF = ReferenceFrame(p->pWheel);
	p->WheelAxle = HP.GetVecRel(RF);
	
	/*
	 * leggo il terreno
	 */
	p->pGround = pDM->ReadNode<StructNode, Node::STRUCTURAL>(HP);
	
	/*
	 * leggo posizione ed orientazione del terreno nel sistema del nodo
	 */
	RF = ReferenceFrame(p->pGround);
	p->GroundPosition = HP.GetPosRel(RF);
	p->GroundDirection = HP.GetVecRel(RF);

	/*
	 * normalizzo l'orientazione del terreno
	 */
	doublereal d = p->GroundDirection.Dot();
	if (d <= std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("Wheellugre(" << pEl->GetLabel() << "): "
			"null direction at line " << HP.GetLineData()
			<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	p->GroundDirection /= sqrt(d);

	/*
	 * Geometria ruota
	 */
	p->dRadius = HP.GetReal();
	
	/*
	 * termine per il calcolo di Delta L
	 */
	p->dRNP = p->dRadius+p->GroundPosition*p->GroundDirection;

	p->dVn = 0.;

	p->F_function = ParseScalarFunction(HP, pDM);
	p->dHystVRef = HP.GetReal();


	/*7
	 * Attrito
	 */
	p->fSlip = false;

	p->fc = 0;

	if (HP.IsKeyWord("slip")) {
		p->fSlip = true;
		
		if (HP.IsKeyWord("preload")) {
			p->Preload = HP.GetReal();
		}
		
		p->fc = ParseFriction(HP,pDM);

		p->kappa0 = HP.GetReal();

		p->re_function = ParseScalarFunction(HP, pDM);

		/*
		 * Parametri di attrito
		 */
		p->pMuY0 = HP.GetDriveCaller();
		p->pMuY1 = HP.GetDriveCaller();
	
		p->dAlphaThreshold = 0.;
		if (HP.IsKeyWord("threshold")) {
			p->dAlphaThreshold = HP.GetReal();
			if (p->dAlphaThreshold < 0.) {
				silent_cerr("Wheellugre(" << pEl->GetLabel() << "): "
					"illegal slip angle threshold " << p->dAlphaThreshold
					<< " at line " << HP.GetLineData() << std::endl);
				p->dAlphaThreshold = fabs(p->dAlphaThreshold);
			}
		}
	}

	p->kb = 0;
	p->pRim = 0;
	if (HP.IsKeyWord("damping")) {
                p->pRim = pDM->ReadNode<StructNode, Node::STRUCTURAL>(HP);
		p->kb = HP.GetReal();
		if (p->kb < 0.) {
			silent_cerr("Wheel2(" << pEl->GetLabel() << "): "
				"illegal tread stiffness " << p->kb
				<< " at line " << HP.GetLineData() << std::endl);
			p->kb = fabs(p->kb);
		}
	}
	
	std::ostream& out = pDM->GetLogFile();
	out << "Wheellugre: " << pEl->GetLabel()
		<< " " << p->pWheel->GetLabel()	//node label
		<< " " << p->WheelAxle		//wheel axle
		<< " " << p->pGround->GetLabel()//ground label
		<< " " << p->GroundDirection	//ground direction
		<< " " << p->dRadius		//wheel radius
		<< std::endl;
	
	return (void *)p;
}

#if 0
static unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_dof");
	return 0;
}
#endif

static unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_dof");
	module_wheellugre* p = (module_wheellugre*)pEl->pGetData();
	if (p->fSlip) {
		return p->fc->iGetNumDof();
	}
	return 0;
}

static DofOrder::Order
set_dof(const LoadableElem*pEl, unsigned int i)
{
	module_wheellugre* p = (module_wheellugre*)pEl->pGetData();
	if (p->fSlip) {
		return p->fc->GetDofType(i);
	}
	return DofOrder::ALGEBRAIC;
}

static DofOrder::Order
set_eq(const LoadableElem*pEl, unsigned int i)
{
	module_wheellugre* p = (module_wheellugre*)pEl->pGetData();
	if (p->fSlip) {
		return p->fc->GetEqType(i);
	}
	return DofOrder::ALGEBRAIC;
}

static void
output(const LoadableElem* pEl, OutputHandler& OH)
{
	DEBUGCOUTFNAME("output");
	
	if (pEl->fToBeOutput()) {
		module_wheellugre* p = (module_wheellugre *)pEl->pGetData();      
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << pEl->GetLabel()	/* 1:	label */
			<< " ", p->F.Write(out, " ")	/* 2-4:	force */
			<< " ", p->M.Write(out, " ")	/* 5-7:	moment */
			<< " " << p->dInstRadius	/* 8:	inst. radius */
			<< " " << p->dDeltaL		/* 9:	radial deformation */
			<< " " << p->dVn 		/* 10:	radial deformation velocity */
			<< " " << p->dSr		/* 11:	slip ratio */
			<< " " << 180./M_PI*p->dAlpha	/* 12:	slip angle */
			<< " " << p->dMuX		/* 13:	longitudinal friction coefficient */
			<< " " << p->dMuY		/* 14:	lateral friction coefficient */
			<< " " << p->dVa		/* 15:	axis relative velocity */
			<< " " << p->dVc		/* 16:	POC relative velocity */
			<< " " << p->V_cx		/* 17:	 */
			<< " " << p->V_cr		/* 18:	 */
			<< " " << p->V_c_sx		/* 19:	 */
			<< " " << p->zeta		/* 20:	 */
			<< " " << p->zetap		/* 21:	 */
			<< " " << p->re			/* 22:	 */
			<< " " << p->fwd		/* 23-25:*/
			<< " " << p->dFn[0]		/* 26:	 */
			<< " " << p->a			/* 27:	 */
			<< " " << p->mu			/* 28:	 */
			<< " " << p->F_cx		/* 29:	 */
			<< " " << p->kappa		/* 30:	 */
			<< " " << p->L;			/* 31:	 */
		if (p->fSlip) {
			  p->v[0] = p->V_c_sx;
		          out << " " << p->fc->fc();	/* 30:	 */
			  // out << p->fc->Output(out,p->dFn,p->v,p->params,p->zeta,p->zetap);
      		}
		out << std::endl;

	}
}

#if 0
static std::ostream&
restart(const LoadableElem* pEl, std::ostream& out)
{
	DEBUGCOUTFNAME("restart");
	return out << "not implemented yet;" << std::endl;
}
#endif

static void
work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
	DEBUGCOUTFNAME("work_space_dim");
	module_wheellugre* p = (module_wheellugre *)pEl->pGetData();
	if (p->fSlip) {
		*piNumRows = 18 + 1;
		*piNumCols = 18 + 1;
	} else {
		*piNumRows = 18;
		*piNumCols = 18;
	}
}

static VariableSubMatrixHandler& 
ass_jac(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
	DEBUGCOUTFNAME("ass_jac");   

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	
	module_wheellugre* p = (module_wheellugre*)pEl->pGetData();

	/*
	 * Dimensiono la matrice
	 */
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);

	WM.ResizeReset(iNumRows, iNumCols);

	integer iGroundFirstPosIndex = p->pGround->iGetFirstPositionIndex();
	integer iGroundFirstMomIndex = p->pGround->iGetFirstMomentumIndex();
	integer iWheelFirstPosIndex = p->pWheel->iGetFirstPositionIndex();
	integer iWheelFirstMomIndex = p->pWheel->iGetFirstMomentumIndex();
	integer iRimFirstPosIndex = iWheelFirstPosIndex;
	integer iRimFirstMomIndex = iWheelFirstMomIndex;
	if (p->pRim) {
		iRimFirstPosIndex = p->pRim->iGetFirstPositionIndex();
		iRimFirstMomIndex = p->pRim->iGetFirstMomentumIndex();
	}

	integer iFirstReactionIndex = pEl->iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iGroundFirstMomIndex+iCnt);
		WM.PutRowIndex(6 + iCnt, iWheelFirstMomIndex + iCnt);
		WM.PutRowIndex(12+iCnt, iRimFirstMomIndex+iCnt);

		WM.PutColIndex(iCnt, iGroundFirstPosIndex+iCnt);
		WM.PutColIndex(6 + iCnt, iWheelFirstPosIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iRimFirstPosIndex + iCnt);
	}
	if (p->fSlip) {
		for (unsigned int i = 1; i <= pEl->iGetNumDof(); i++) {
			WM.PutRowIndex(18 + i, iFirstReactionIndex + i);
			WM.PutColIndex(18 + i, iFirstReactionIndex + i);
		}
	
		if (p->dDeltaL <= 0.) {
			WM.PutCoef(18 + 1, 18 + 1, -dCoef);
		} else {

			ExpandableRowVector dfc;
			ExpandableRowVector dF;
			ExpandableRowVector dv;
			//variation of reaction force
			dF.ReDim(0);
			//variation of relative velocity
			dv.ReDim(0);
			p->v[0] = p->V_c_sx;
			// p->fc->AssJac(WM,dfc,18,iFirstReactionIndex,dCoef,p->dFn,p->v,p->params,XCurr,XPrimeCurr,dF,dv);
			p->fc->AssJac(WM,dfc,18,iFirstReactionIndex,dCoef,p->dFn[0],p->v[0],XCurr,XPrimeCurr,dF,dv);
		}
	}

	
	return WorkMat;
}

static SubVectorHandler& 
ass_res(LoadableElem* pEl, 
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ass_res");
	
	module_wheellugre* p = (module_wheellugre *)pEl->pGetData();

	/*
	 * Orientazione del terreno nel sistema assoluto
	 */
	Vec3 n = p->pGround->GetRCurr()*p->GroundDirection;
	
	/*
	 * Distanza Wheel Ground nel sistema assoluto
	 */
	Vec3 x = p->pWheel->GetXCurr()-p->pGround->GetXCurr();

	/*
	 * Smorzamento kb
	 */
	Vec3 xWheel_xRim(0., 0., 0.);
	if (p->pRim) {
		xWheel_xRim = p->pRim->GetXCurr()-p->pGround->GetXCurr() - x;
	}

	/*
	 * se dDeltaL > 0 c'e' contatto, altrimenti no
	 */
	p->dDeltaL = p->dRNP - x*n;

	/*
	 * Reset dati per output
	 */
	p->dInstRadius = p->dRadius-p->dDeltaL;
	
	p->dSr = 0.;
	p->dAlpha = 0.;

	p->dMuX = 0.;
	p->dMuY = 0.;

	p->dVa = 0.;
	p->dVc = 0.;

	integer iFirstReactionIndex = pEl->iGetFirstIndex();
	doublereal zeta = XCurr.dGetCoef(iFirstReactionIndex + 1);
	doublereal zetap = XPrimeCurr.dGetCoef(iFirstReactionIndex + 1);

	p->zeta = zeta;
	p->zetap = zetap;

	if (p->dDeltaL <= 0.) {
		
		p->F = Zero3;
		p->M = Zero3;

		p->dInstRadius = p->dRadius;
		p->dDeltaL = 0.;
		
		/*
		 * Non assemblo neppure il vettore ;)
		 */
		WorkVec.Resize(0);

		if (p->fSlip) {
			WorkVec.ResizeReset(pEl->iGetNumDof());
			for (unsigned int i = 1; i <= pEl->iGetNumDof(); i++) {
				WorkVec.PutRowIndex(i, iFirstReactionIndex + i);
			}
			WorkVec.IncCoef(1, zeta);
		}
		
		return WorkVec;
	}
	
	/*
	 * Dimensiono il vettore
	 */
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
	WorkVec.ResizeReset(iNumRows);

	integer iGroundFirstMomIndex = p->pGround->iGetFirstMomentumIndex();
	integer iWheelFirstMomIndex = p->pWheel->iGetFirstMomentumIndex();
	integer iRimFirstMomIndex = iWheelFirstMomIndex;
	if (p->pRim) {
		iRimFirstMomIndex = p->pRim->iGetFirstMomentumIndex();
	}

	/*
	 * Indici equazioni
	 */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iGroundFirstMomIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iWheelFirstMomIndex + iCnt);
		WorkVec.PutRowIndex(12+iCnt, iRimFirstMomIndex+iCnt);
	}
	if (p->fSlip) {
		for (unsigned int i = 1; i <= pEl->iGetNumDof(); i++) {
			WorkVec.PutRowIndex(18 + i, iFirstReactionIndex + i);
		}
	}


	/*
	 * Velocita' tra Wheel (nell'asse)
	 * e Ground nel sistema assoluto
	 */
	Vec3 va = p->pWheel->GetVCurr()
		-p->pGround->GetVCurr()-(p->pGround->GetWCurr()).Cross(
			p->pGround->GetRCurr()*p->GroundPosition
			);

	p->dVa = (va - n*(n*va)).Norm();
	
	/*
	 * Velocita' tra Wheel (nel punto di contatto) 
	 * e Ground nel sistema assoluto
	 */
	Vec3 v = va - (p->pWheel->GetWCurr()).Cross(n*p->dInstRadius);
	
	/*
	 * Componente normale al terreno della velocita'
	 * (positiva se la ruota si allontana dal terreno)
	 */
	p->dVn = n*v;
	
	p->dVc = (v - n*p->dVn).Norm();
	
	/*
	 * Il punto di applicazione della forza e' xW - R * nG ;
	 * per la forza normale si puo' usare anche la posizione
	 * della ruota, nell'ipotesi che il punto di contatto
	 * stia nell'intersezione della normale al ground passante
	 * per l'asse ruota.
	 * 
	 * pc e' corretto per il raggio effettivo (diverso da dInstRadius)
	 * se c'e' attrito.
	 *
	 * FIXME: perche' dRadius invece di dInstRadius?
	 */
	Vec3 pc = p->pWheel->GetXCurr() - (n*p->dInstRadius);

	/*
	 * Forza verticale
	 */
	p->dFn[0] = (*p->F_function)(p->dDeltaL);
	p->F = n * (p->dFn[0]* (1. - tanh(p->dVn/p->dHystVRef)));
		

	p->sigma_c = 0.;
	p->V_c_sx = 0.;
	p->V_cr = 0.;
	
	bool ChangeJac(false);
			
	if (p->fSlip) {
		/*
		 * compute contact length: ESDU 86005
		 *
		*/
		{
			doublereal d = 2. * p->dRadius;
			p->L = d * 1.7 * std::sqrt(p->dDeltaL / d - std::pow(p->dDeltaL / d, 2));
			// parachute for slight contact
			p->L = std::max(0.05 * p->dRadius, p->L);
			p->kappa = p->kappa0 / p->L;
		}
	
		p->dFn[0] = std::max(p->dFn[0], p->Preload);
		/*
		 * Direzione di "avanzamento": asse ruota cross normale
		 * al ground
		 */
		p->fwd = (p->pWheel->GetRCurr()*p->WheelAxle).Cross(n);
		doublereal d = p->fwd.Dot();
		if (d < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("Wheellugre(" << pEl->GetLabel() << "): "
				"wheel axle is (nearly) orthogonal "
				"to the ground" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		p->fwd /= sqrt(d);


		/*
		 * Slip velocity
		 */
		p->re = (*p->re_function)(p->dFn[0]);
		/* correggo pc per il raggio effettivo !*/
		pc = p->pWheel->GetXCurr() - (n*p->re);
		
		p->V_cr = p->fwd.Dot((p->pWheel->GetWCurr()).Cross(n*p->re));
		p->V_cx = p->fwd.Dot(va);
		
		p->V_c_sx = p->V_cx - p->V_cr;
		
		p->dSr = std::abs(p->V_c_sx)/(std::abs(p->V_cx)+0.000001);
// 		if (p->dSr > 1.) {
// 			p->dSr = 1.;
// 		}
		
		/* calcolo F_cx */
		//FIXME: TODO
		p->v[0] = p->V_c_sx;
		p->v[1] = p->V_cr;
// 		p->params[0] = p->kappa/p->re;
		p->params[0] = p->kappa;
		try {
			// p->fc->AssRes(WorkVec, 18, iFirstReactionIndex, p->dFn, p->v, p->params, XCurr, XPrimeCurr);
			p->fc->AssRes(WorkVec, 18, iFirstReactionIndex, p->dFn[0], p->v[0], XCurr, XPrimeCurr);
		}
		catch (Elem::ChangedEquationStructure&) {
			ChangeJac = true;
		}
		p->mu = p->fc->fc();
	
		/*
		 * Direzione laterale: normale cross forward
		 */
		Vec3 lat = n.Cross(p->fwd);

		/*
		 * Velocita' laterale del mozzo
		 */
		doublereal dvay = lat.Dot(va);

		/*
		 * Angolo di deriva del mozzo
		 * Nota: ristretto al Io-IVo quadrante
		 * questa threshold sul modulo della velocita' fa in modo
		 * che l'angolo vada a zero se il modulo della velocita'
		 * e' troppo piccolo
		 */
		p->dAlpha = atan2(dvay, fabs(p->V_cx));
		if (p->dAlphaThreshold > 0.) {
			doublereal dtmp = tanh(sqrt(p->V_cx*p->V_cx + dvay*dvay)/p->dAlphaThreshold);
			p->dAlpha *= dtmp * dtmp;
		}

		/*
		 * Coefficiente di attrito longitudinale
		 */
		/*
		 * Nota: alpha/(pi/2) e' compreso tra -1. e 1.
		 */
		p->dMuX = p->mu*(1. - fabs(p->dAlpha)/M_PI_2);
		
		/*
		 * Correggo le forze:
		 * uso il coefficiente di attrito longitudinale
		 * con il segno della velocita' del punto di contatto
		 */
		//if (std::abs(p->V_c_sx) > p->dvThreshold) {
		p->F -= p->fwd*p->dMuX*p->dFn[0];
				//*copysign(1.,p->V_c_sx);
		//} else {
		//	p->F += p->fwd*p->F_cx*(1. - fabs(p->dAlpha)/M_PI_2)
		//		*copysign(1.,p->zeta);
		//}

		if (dvay != 0.) {
			doublereal dAlpha = p->dAlpha;

			doublereal dMuY0 = p->pMuY0->dGet(dAlpha);
			doublereal dMuY1 = p->pMuY1->dGet(dAlpha);
			
			p->dMuY = dMuY0 + (dMuY1 - dMuY0)*p->dSr;

			/*
			 * Correggo le forze
			 */
			p->F -= lat*p->dFn[0]*p->dMuY;
		}
	}

	/*
	 * Momento
	 */
	p->M = (pc - p->pWheel->GetXCurr()).Cross(p->F);

	/*
	 * Smorzamento
	 */
	p->Fdamp = p->pWheel->GetWCurr().Cross(xWheel_xRim)*p->kb;

	WorkVec.Sub(1, p->F);
	WorkVec.Sub(4, (pc - p->pGround->GetXCurr()).Cross(p->F));
	WorkVec.Add(7, p->F);
	WorkVec.Add(10, p->M);

	WorkVec.Add(7, p->Fdamp);
	WorkVec.Sub(13, p->Fdamp);
	
	if (ChangeJac) {
		throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
	}

	return WorkVec;
}

#if 0
static void
before_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev)
{
   DEBUGCOUTFNAME("before_predict");
}

static void
after_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP)
{
	DEBUGCOUTFNAME("after_predict");
}

static void
update(LoadableElem* pEl, 
		const VectorHandler& X,
		const VectorHandler& XP)
{
	DEBUGCOUTFNAME("update");
}

static unsigned int
i_get_initial_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_initial_num_dof");
	return 0;
}

static void
initial_work_space_dim(const LoadableElem* pEl, 
		integer* piNumRows, 
		integer* piNumCols)
{
	DEBUGCOUTFNAME("initial_work_space_dim");
	*piNumRows = 0;
	*piNumCols = 0;
}

static VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("initial_ass_jac");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(iNumRows, iNumCols);
#if 0	
	module_wheellugre* p = (module_wheellugre *)pEl->pGetData();
#endif /* 0 */
   
	/* set sub-matrix indices and coefs */

	return WorkMat;
}

static SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("initial_ass_res");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
	
	WorkVec.Resize(iNumRows);

#if 0
	module_wheellugre* p = (module_wheellugre *)pEl->pGetData(); 
#endif /* 0 */
	
	/* set sub-vector indices and coefs */
   
	return WorkVec;
}

static void
set_value(const LoadableElem* pEl, DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	DEBUGCOUTFNAME("set_value");
}
   
static void
set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
	DEBUGCOUTFNAME("set_initial_value");
}
#endif

static unsigned int
i_get_num_priv_data(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_priv_data");
	return 0;
}

static void
destroy(LoadableElem* pEl)
{
	DEBUGCOUTFNAME("destroy");
	module_wheellugre* p = (module_wheellugre *)pEl->pGetData();
	
	SAFEDELETE(p);
}

static int
i_get_num_connected_nodes(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_connected_nodes");
	module_wheellugre* p = (module_wheellugre *)pEl->pGetData();
	
	/* wheel + ground */
	if (p->pRim) {
		return 3;
	}
	return 2;
}

static void
get_connected_nodes(const LoadableElem* pEl, 
		std::vector<const Node *>& connectedNodes)
{
	DEBUGCOUTFNAME("get_connected_nodes");
	module_wheellugre* p = (module_wheellugre *)pEl->pGetData();
	
	/* wheel + ground */
	if (p->pRim) {
		connectedNodes.resize(3);
		connectedNodes[2] = p->pRim;
	} else {
		connectedNodes.resize(2);
	}

	connectedNodes[0] = p->pWheel;
	connectedNodes[1] = p->pGround;
}

#ifdef MBDYN_MODULE
static
#endif /* MBDYN_MODULE */
struct LoadableCalls module_wheellugre_lc = {
    LOADABLE_VERSION_SET(1, 7, 0), // loadable_version

    "wheellugre",                   // name
    "1.0.0",                        // version
    "Dipartimento di Ingegneria Aerospaziale, Politecnico di Milano", // vendor
    "tire force model for landing gear analysis\n", // description

    read,                           // read

    i_get_num_dof,                  // i_get_num_dof
    set_dof,                        // set_dof
    set_eq,                         // set_eq
    output,                         // output
    NULL,                           // restart
    work_space_dim,                 // work_space_dim
    ass_jac,                        // ass_jac
    NULL,                           // ass_mats
    ass_res,                        // ass_res
    NULL,                           // before_predict
    NULL,                           // after_predict
    NULL,                           // update
    NULL,                           // after_convergence
    NULL,                           // i_get_initial_num_dof
    NULL,                           // initial_work_space_dim
    NULL,                           // initial_ass_jac
    NULL,                           // initial_ass_res
    NULL,                           // set_value
    NULL,                           // set_initial_value
    i_get_num_priv_data,            // i_get_num_priv_data
    NULL,                           // i_get_priv_data_idx
    NULL,                           // d_get_priv_data
    i_get_num_connected_nodes,      // i_get_num_connected_nodes
    get_connected_nodes,            // get_connected_nodes

    destroy                         // destroy
};
extern "C" {

//void *wheellugre = &module_wheellugre_lc;
void *calls = &module_wheellugre_lc;

#ifndef STATIC_MODULES
extern "C" int
module_init(const char *s, void *dm, void *)
{
	DataManager *pDM = (DataManager *)dm;
	
	pDM->SetLoadableElemModule(module_wheellugre_lc.name, &module_wheellugre_lc);

	return 0;
}
#endif /* !STATIC_MODULES */

}


