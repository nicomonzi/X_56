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

/* Giunti sferici */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "spherj.h"
#include "Rot.hh"
#include "SphericalQR.h"


/* SphericalHingeJoint - begin */

const unsigned int SphericalHingeJoint::NumSelfDof(3);
const unsigned int SphericalHingeJoint::NumDof(15);

/* Costruttore non banale */
SphericalHingeJoint::SphericalHingeJoint(unsigned int uL, const DofOwner* pDO,
					 const StructNode* pN1, 
					 const StructNode* pN2,
					 const Vec3& dTmp1, const Mat3x3& RTmp1h,
					 const Vec3& dTmp2, const Mat3x3& RTmp2h,
					 const OrientationDescription& od_a,
					 flag fOut,
                     const doublereal rr,
                     const doublereal pref,
                     BasicShapeCoefficient2D *const sh,
                     BasicFriction2D *const f)
: Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2), 
d1(dTmp1), R1h(RTmp1h),
d2(dTmp2), R2h(RTmp2h), 
F(Zero3),
Sh_c(sh), fc(f), preF(pref), r(rr), Ffrict1(Zero3), Ffrict2(Zero3), M1(Zero3), M2(Zero3),
Q(Eye3), Qold(Eye3), Fold(Zero3), reset_Q(true), compute_Q(false),
od(od_a)
{
   NO_OP;
}


/* Distruttore banale */
SphericalHingeJoint::~SphericalHingeJoint(void)
{
	if (Sh_c) {
		delete Sh_c;
	}

	if (fc) {
		delete fc;
	}
};

std::ostream&
SphericalHingeJoint::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"reaction forces [Fx,Fy,Fz]" << std::endl;

	if (bInitial) {
		iIndex += NumSelfDof;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"reaction force derivatives [FPx,FPy,FPz]" << std::endl;
	}

	iIndex += NumSelfDof;
	if (fc) {
		integer iFCDofs = fc->iGetNumDof();
		if (iFCDofs > 0) {
			out << prefix << iIndex + 1;
			if (iFCDofs > 1) {
				out << "->" << iIndex + iFCDofs;
			}
			out << ": friction dof(s)" << std::endl
				<< "        ", fc->DescribeDof(out, prefix, bInitial);
		}
	}

	return out;
}

static const char xyz[] = "xyz";

void
SphericalHingeJoint::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "SphericalHingeJoint(" << GetLabel() << ")";

	unsigned short nself = NumSelfDof;
	if (bInitial) {
		nself *= 2;
	}
	desc.resize(nself);
	if (fc && (i == -1 || i >= nself)) {
		fc->DescribeDof(desc, bInitial, i - nself);
		if (i != -1) {
			desc[0] = os.str() + ": " + desc[0];
			desc[1] = os.str() + ": " + desc[1];
			return;
		}
	}

	if (i == -1) {
		// move fc desc to the end
		unsigned short nfc = 0;
		if (fc) {
			nfc = desc.size();
		}
		desc.resize(nfc + nself);
		for (unsigned ii = nfc; ii-- > 0; ) {
			desc[nself + ii] = os.str() + ": " + desc[ii];
		}

		std::string name = os.str();

		for (unsigned ii = 0; ii < 3; ii++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": reaction force f" << xyz[ii];
			desc[ii] = os.str();
		}

		if (bInitial) {
			for (unsigned ii = 0; ii < 3; ii++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": reaction force derivative fP" << xyz[ii];
				desc[3 + ii] = os.str();
			}
		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
		case 2:
			os << ": reaction force f" << xyz[i];
			break;

		case 3:
		case 4:
		case 5:
			os << ": reaction force derivative fP" << xyz[i - 3];
			break;
		}
		desc[0] = os.str();
	}
}

std::ostream&
SphericalHingeJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"position constraints [Px1=Px2,Py1=Py2,Pz1=Pz2]" << std::endl;

	if (bInitial) {
		iIndex += NumSelfDof;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"velocity constraints [vx1=vx2,vy1=vy2,vz1=vz2]" << std::endl;
	}

	iIndex += NumSelfDof;
	if (fc) {
		integer iFCDofs = fc->iGetNumDof();
		if (iFCDofs > 0) {
			out << prefix << iIndex + 1;
			if (iFCDofs > 1) {
				out << "->" << iIndex + iFCDofs;
			}
			out << ": friction equation(s)" << std::endl
				<< "        ", fc->DescribeEq(out, prefix, bInitial);
		}
	}

	return out;
}

void
SphericalHingeJoint::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "SphericalHingeJoint(" << GetLabel() << ")";

	unsigned short nself = NumSelfDof;
	if (bInitial) {
		nself *= 2;
	}
	if (fc && (i == -1 || i >= nself)) {
		fc->DescribeEq(desc, bInitial, i - nself);
		if (i != -1) {
			desc[0] = os.str() + ": " + desc[0];
			return;
		}
	}

	if (i == -1) {
		// move fc desc to the end
		unsigned short nfc = 0;
		if (fc) {
			nfc = desc.size();
		}
		desc.resize(nfc + nself);
		for (unsigned ii = nfc; ii-- > 0; ) {
			desc[nself + ii] = os.str() + ": " + desc[ii];
		}

		std::string name = os.str();

		for (unsigned ii = 0; ii < 3; ii++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": position constraint P" << xyz[ii];
			desc[ii] = os.str();
		}

		if (bInitial) {
			for (unsigned ii = 0; ii < 3; ii++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": position constraint derivative v" << xyz[ii];
				desc[3 + ii] = os.str();
			}

		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
		case 2:
			os << ": position constraint P" << xyz[i];
			break;

		case 3:
		case 4:
		case 5:
			os << ": position constraint derivative v" << xyz[i - 3];
			break;

		}
		desc[0] = os.str();
	}
}

void SphericalHingeJoint::AfterConvergence(const VectorHandler& X,
				const VectorHandler& XP) {
	if (fc) {
		//FIXME
		//doublereal v = (Omega1-Omega2).Dot(e3a)*r;
		Qold = Q;
		Fold = F;
		Vec3 Omega1(pNode1->GetWCurr());
		Vec3 Omega2(pNode2->GetWCurr());
		Vec3 Omegar = Omega1 - Omega2;
		F = Vec3(X, iGetFirstIndex()+1);
		doublereal modF = F.Norm();
		if (modF <= preF / 2.) {
			// std::cout << "Condition 3" << std::endl;
			compute_Q = false;
			// std::cout << "reset_Q: " << reset_Q << "; compute_Q: " << compute_Q << std::endl;
		} else {
			SphericalQR(F, Q, !(reset_Q), Qold);
		}
		d2D v;
		v.x[0] = (-Q.GetCol(1)).Cross(Omegar).Dot(Q.GetCol(2))*r;
		v.x[1] = (-Q.GetCol(1)).Cross(Omegar).Dot(Q.GetCol(3))*r;
		//reaction norm
		// doublereal modF = std::max(modF, preF);
		fc->AfterConvergence(modF, v, X, XP, iGetFirstIndex()+NumSelfDof);
	}
}

/* Contributo al file di restart */
std::ostream& SphericalHingeJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", spherical hinge, "
     << pNode1->GetLabel() << ", reference, node, ",
     d1.Write(out, ", ")  << ", hinge, reference, node, 1, ", (R1h.GetCol(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetCol(2)).Write(out, ", ") << ", "
     << pNode2->GetLabel() << ", reference, node, ",
     d2.Write(out, ", ") << ", hinge, reference, node, 1, ", (R2h.GetCol(1)).Write(out, ", ")
     << ", 2, ", (R2h.GetCol(2)).Write(out, ", ") << ';' << std::endl;
   
   return out;
}

void SphericalHingeJoint::Restart(RestartData& oData, RestartData::RestartAction eAction)
{
     oData.Sync(RestartData::ELEM_JOINTS, GetLabel(), "d1", d1, eAction);
     oData.Sync(RestartData::ELEM_JOINTS, GetLabel(), "d2", d2, eAction);
     oData.Sync(RestartData::ELEM_JOINTS, GetLabel(), "F", F, eAction);
}

/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
SphericalHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			    doublereal dCoef,
			    const VectorHandler& XCurr ,
			    const VectorHandler& XPrimeCurr )
{
   DEBUGCOUT("Entering SphericalHingeJoint::AssJac()" << std::endl);
      
   /* Setta la sottomatrice come piena (e' un po' dispersivo, ma lo jacobiano
    * e' complicato */
   FullSubMatrixHandler& WM = WorkMat.SetFull();

   /* Ridimensiona la sottomatrice in base alle esigenze */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici delle varie incognite */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

      /* Setta gli indici delle equazioni */
   for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }

   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }

   Vec3 dTmp1(pNode1->GetRRef()*d1);
   Vec3 dTmp2(pNode2->GetRRef()*d2);
   
   
   /* 
    * L'equazione di vincolo afferma che il punto in cui si trova la
    * cerniera deve essere consistente con la posizione dei due nodi:
    *      x2 + d2 = x1 + d1
    * 
    * con: d2 = R2 * d2_0
    *      d1 = R1 * d1_0
    * 
    * La forza e' data dalla reazione vincolare F, nel sistema globale
    * La coppia dovuta all'eccentricita' e' data rispettivamente da:
    *     -d1 /\ F    per il nodo 1,
    *      d2 /\ F    per il nodo 2
    *
    * 
    *         x1   g1        x2     g2        F
    * Q1 |  0      0         0      0         I    | | x1 |   | -F           |
    * G1 |  0      cF/\d1/\  0      0         d1/\ | | g1 |   | -d1/\F       |
    * Q2 |  0      0         0     -cF/\d2/\ -I    | | x2 | = |  F           |
    * G2 |  0      0         0      0        -d2/\ | | g2 |   |  d2/\F       |
    * F  | -c*I    c*d1/\    c*I   -c*d2/\    0    | | F  |   |  x1+d1-x2-d2 |
    * 
    * con d1 = R1*d01, d2 = R2*d02, c = dCoef
    */

   /* Moltiplico la forza per il coefficiente del metodo.
    * Nota: F, la reazione vincolare, e' stata aggiornata da AssRes */
   Vec3 FTmp = F*dCoef;
   
   /* termini di reazione sul nodo 1 */
   WM.PutDiag(1, 13, 1.);
   WM.PutCross(4, 13, dTmp1);
   
   WM.Put(4, 4, Mat3x3(MatCrossCross, FTmp, dTmp1));

   /* termini di reazione sul nodo 2 */
   WM.PutDiag(7, 13, -1.);
   WM.PutCross(10, 13, -dTmp2);

   WM.Put(10, 10, Mat3x3(MatCrossCross, FTmp, -dTmp2));
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* termini di vincolo dovuti al nodo 1 */
   WM.PutDiag(13, 1, -1.);
   WM.PutCross(13, 4, dTmp1);
      
   /* termini di vincolo dovuti al nodo 1 */
   WM.PutDiag(13, 7, 1.);
   WM.PutCross(13, 10, -dTmp2);
   if (fc) {
      //retrive
          //friction coef
      d2D f = fc->fc();
          //shape function
      d2D shc = Sh_c->Sh_c();
          //omega and omega rif
      const Vec3& Omega1(pNode1->GetWCurr());
      const Vec3& Omega2(pNode2->GetWCurr());
      const Vec3 Omegar(Omega1-Omega2);
	  const Vec3& Omega1Ref(pNode1->GetWRef());
	  const Vec3& Omega2Ref(pNode2->GetWRef());
      d2D v;
      v.x[0] = (-Q.GetCol(1)).Cross(Omegar).Dot(Q.GetCol(2))*r;
      v.x[1] = (-Q.GetCol(1)).Cross(Omegar).Dot(Q.GetCol(3))*r;
      // v.x[0] = (Omegar).Dot(Q.GetCol(2))*r;
      // v.x[1] = (Omegar).Dot(Q.GetCol(3))*r;
      ExpandableMatrix dF, dshc, dfc, dv, dQ1, dQ2, dQ3, dOmegar;
      ExpandableRowVector dmodF;
      ExpandableRowVector drealmodF;

      // dQ1, dmodF
      doublereal modF = std::max(F.Norm(), preF);
      dmodF.ReDim(3);
      drealmodF.ReDim(3);
      dQ1.ReDim(3, 2);
      dQ1.SetBlockDim(1, 3);
      dQ1.SetBlockIdx(1, 12+1);
      dQ1.SetBlockDim(2, 1);
    //   if (F.Norm() == 0.) {
		  // // std::cout << "!compute_Q\n";
    //       dmodF.Set(Zero3, 1, 12+1);
    //       dQ1.Set(Zero3x3, 1, 1, 1);
    //       dQ1.SetCol(Zero3, 1, 2, 1);
    //       dQ1.Link(2, &dmodF);
    //   } else {
		  // std::cout << "compute_Q\n";
	  if (F.Norm() < preF){
		  dmodF.Set(Zero3, 1, 12+1);
	  } else {
          dmodF.Set(F/modF, 1, 12+1);
      }
      if (compute_Q) {
          // std::cout << "HERE 1 --------------------" << std::endl;
          doublereal real_modF = F.Norm();
		  // std::cerr << F.Norm() << " " << GetLabel() << std::endl;
          drealmodF.Set(F/real_modF, 1, 12+1);
          dQ1.Set(Eye3/real_modF, 1, 1, 1);
          dQ1.SetCol(-F/(real_modF*real_modF), 1, 2, 1);
          dQ1.Link(2, &drealmodF);
      } else {
          // std::cout << "THERE ********************** " << GetLabel() << std::endl;
		  drealmodF.Set(Zero3, 1, 12+1);
          dQ1.Set(Zero3x3, 1, 1, 1);
          dQ1.SetCol(Zero3, 1, 2, 1);
          dQ1.Link(2, &drealmodF);
      }

      // dQ2
      dQ2.ReDim(3, 1);
      dQ2.SetBlockDim(1, 3);
      dQ2.Set(-Q.GetCol(1).Tens(Q.GetCol(2)), 1, 1, 1);
      dQ2.Link(1, &dQ1);
      // dQ3
      dQ3.ReDim(3, 1);
      dQ3.SetBlockDim(1, 3);
      dQ3.Set(-Q.GetCol(1).Tens(Q.GetCol(3)), 1, 1, 1);
      dQ3.Link(1, &dQ1);

      // dOmegar
	  // Vec3 Omegar(Omega1 - Omega2);
      dOmegar.ReDim(3, 2);
      dOmegar.SetBlockDim(1, 3); // dOmegar/dOmega1
      dOmegar.SetBlockIdx(1, 4);
      dOmegar.Set(Eye3 - Mat3x3(MatCross, Omega1Ref*dCoef), 1, 1, 1);
      dOmegar.SetBlockDim(2, 3); // dOmegar/dOmega2
      dOmegar.SetBlockIdx(2, 10);
      dOmegar.Set(-Eye3 + Mat3x3(MatCross, Omega2Ref*dCoef), 1, 2, 1);

	  // v.x[0] = (-Q.GetCol(1)).Cross(Omegar).Dot(Q.GetCol(2))*r;
	  // v.x[1] = (-Q.GetCol(1)).Cross(Omegar).Dot(Q.GetCol(3))*r;
	  // dv
      dv.ReDim(2, 4);
      dv.SetBlockDim(1, 3); // dv/dOmegar
      dv.SetRow(Q.GetCol(1).Cross(Q.GetCol(2))*r, 1, 1, 1);
      dv.SetRow(Q.GetCol(1).Cross(Q.GetCol(3))*r, 2, 1, 1);
      dv.Link(1, &dOmegar);
      dv.SetBlockDim(2, 3); // dv/dQ2
      dv.SetRow(-Q.GetCol(1).Cross(Omegar)*r, 1, 2, 1);
      dv.SetRow(Zero3, 2, 2, 1);
      dv.Link(2, &dQ2);
      dv.SetBlockDim(3, 3); // dv/dQ3
      dv.SetRow(Zero3   , 1, 3, 1);
      dv.SetRow(-Q.GetCol(1).Cross(Omegar)*r, 2, 3, 1);
      dv.Link(3, &dQ3);
	  dv.SetBlockDim(4, 3); //dv/dQ1
	  dv.SetRow(-Omegar.Cross(Q.GetCol(2))*r, 1, 4, 1);
	  dv.SetRow(-Omegar.Cross(Q.GetCol(3))*r, 2, 4, 1);
	  dv.Link(4, &dQ1);

      // dshc
      fc->AssJac(WM, dfc, 12+NumSelfDof, iFirstReactionIndex+NumSelfDof, dCoef, modF, v,
      		XCurr, XPrimeCurr, dmodF, dv);
      Sh_c->dSh_c(dshc, f, modF, v, dfc, dmodF, dv);
      // if (compute_Q) {
		  // std::cout << "Dentro Jac force" << std::endl;
          ExpandableMatrix dF1, dF2;

		  // FIXME: F.Norm() or modF?
		  // Ffrict1 = -Q.GetCol(2)*F.Norm()*f.x[0];
		  // Ffrict2 = -Q.GetCol(3)*F.Norm()*f.x[1];

		  dF1.ReDim(3, 3);
          dF1.SetBlockDim(1, 3); // dF1/dQ2
          dF1.Set(-Eye3 * modF * f.x[0], 1, 1, 1);
          dF1.Link(1, &dQ2);
          dF1.SetBlockDim(2, 1); // dF1/dnormF
          dF1.SetCol(-Q.GetCol(2) * f.x[0], 1, 2, 1);
          dF1.Link(2, &dmodF);
          dF1.SetBlockDim(3, 2); // dF1/dfc
          dF1.SetCol(-Q.GetCol(2) * modF, 1, 3, 1);
          dF1.SetCol(Zero3, 1, 3, 2);
          dF1.Link(3, &dfc);

          dF2.ReDim(3, 3);
          dF2.SetBlockDim(1, 3); // dF2/dQ3
          dF2.Set(-Eye3 * modF * f.x[1], 1, 1, 1);
          dF2.Link(1, &dQ3);
          dF2.SetBlockDim(2, 1); // dF2/dnormF
          dF2.SetCol(-Q.GetCol(3) * f.x[1], 1, 2, 1);
          dF2.Link(2, &dmodF);
          dF2.SetBlockDim(3, 2); // dF1/dfc
          dF2.SetCol(Zero3, 1, 3, 1);
          dF2.SetCol(-Q.GetCol(3) * modF, 1, 3, 2);
          dF2.Link(3, &dfc);

		  // WorkVec.Sub(1, Ffrict1);
		  // WorkVec.Sub(1, Ffrict2);
		  // WorkVec.Add(7, Ffrict1);
		  // WorkVec.Add(7, Ffrict2);
		  dF1.AddTo(WM, 1, 1.);
		  dF2.AddTo(WM, 1, 1.);
		  dF1.SubFrom(WM, 7, 1.);
		  dF2.SubFrom(WM, 7, 1.);

		  // WorkVec.Sub(4, dTmp1.Cross(Ffrict1)); /* Sfrutto  F/\d = -d/\F */
		  // WorkVec.Sub(4, dTmp1.Cross(Ffrict2)); /* Sfrutto  F/\d = -d/\F */
		  // WorkVec.Add(10, dTmp2.Cross(Ffrict1));
		  // WorkVec.Add(10, dTmp2.Cross(Ffrict2));

		  ExpandableMatrix dMF1;
          dMF1.ReDim(3, 1);
          dMF1.SetBlockDim(1, 3);
          dMF1.Set(Mat3x3(MatCross, dTmp1), 1, 1, 1);
          dMF1.Link(1, &dF1);
          dMF1.AddTo(WM, 4);
          dMF1.Link(1, &dF2);
          dMF1.AddTo(WM, 4);

		  ExpandableMatrix dMF2;
          dMF2.ReDim(3, 1);
          dMF2.SetBlockDim(1, 3);
          dMF2.Set(Mat3x3(MatCross, dTmp2), 1, 1, 1);
          dMF2.Link(1, &dF1);
          dMF2.SubFrom(WM, 10);
          dMF2.Link(1, &dF2);
          dMF2.SubFrom(WM, 10);


		  // WM.Put(4, 4, Mat3x3(MatCrossCross, FTmp, dTmp1));
		  // WM.Put(10, 10, Mat3x3(MatCrossCross, FTmp, -dTmp2));
		  WM.Add(4, 4, Mat3x3(MatCrossCross, Ffrict1*dCoef, dTmp1));
		  WM.Add(4, 4, Mat3x3(MatCrossCross, Ffrict2*dCoef, dTmp1));
		  WM.Add(10, 10, Mat3x3(MatCrossCross, Ffrict1*dCoef, -dTmp2));
		  WM.Add(10, 10, Mat3x3(MatCrossCross, Ffrict2*dCoef, -dTmp2));

      // }
      // Vec3 M1 = Q.GetCol(1).Cross(Q.GetCol(2)) * shc.x[0] * r * modF;
      // Vec3 M2 = Q.GetCol(1).Cross(Q.GetCol(3)) * shc.x[1] * r * modF;
      ExpandableMatrix dM1, dM2;

      dM1.ReDim(3, 4);
	  dM1.SetBlockDim(1, 3); //dM1/dQ1
	  dM1.Set(-Mat3x3(MatCross, Q.GetCol(2)) * shc.x[0] * r * modF, 1, 1, 1);
	  dM1.Link(1, &dQ1);

	  dM1.SetBlockDim(2, 3); //dM1/dQ2
	  dM1.Set(Mat3x3(MatCross, Q.GetCol(1)) * shc.x[0] * r * modF, 1, 2, 1);
	  dM1.Link(2, &dQ2);

	  dM1.SetBlockDim(3, 2); //dM1/dshc
	  dM1.SetCol(Q.GetCol(1).Cross(Q.GetCol(2)) * r * modF, 1, 3, 1);
	  dM1.SetCol(Zero3, 1, 3, 2);
	  dM1.Link(3, &dshc);

	  dM1.SetBlockDim(4, 1); //dM1/dmodF
	  dM1.SetCol(Q.GetCol(1).Cross(Q.GetCol(2)) * r * shc.x[0], 1, 4, 1);
	  dM1.Link(4, &dmodF);


      dM2.ReDim(3, 4);
	  dM2.SetBlockDim(1, 3); //dM2/dQ1
	  dM2.Set(-Mat3x3(MatCross, Q.GetCol(3)) * shc.x[1] * r * modF, 1, 1, 1);
	  dM2.Link(1, &dQ1);

	  dM2.SetBlockDim(2, 3); //dM2/dQ3
	  dM2.Set(Mat3x3(MatCross, Q.GetCol(1)) * shc.x[1] * r * modF, 1, 2, 1);
	  dM2.Link(2, &dQ3);

	  dM2.SetBlockDim(3, 2); //dM2/dshc
	  dM2.SetCol(Zero3, 1, 3, 1);
	  dM2.SetCol(Q.GetCol(1).Cross(Q.GetCol(3)) * r * modF, 1, 3, 2);
	  dM2.Link(3, &dshc);

	  dM2.SetBlockDim(4, 1); //dM2/dmodF
	  dM2.SetCol(Q.GetCol(1).Cross(Q.GetCol(3)) * r * shc.x[1], 1, 4, 1);
	  dM2.Link(4, &dmodF);

	  dM1.AddTo(WM, 4, 1.);
	  dM2.AddTo(WM, 4, 1.);
	  dM1.SubFrom(WM, 10, 1.);
	  dM2.SubFrom(WM, 10, 1.);
   }

   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& SphericalHingeJoint::AssRes(SubVectorHandler& WorkVec,
					      doublereal dCoef,
					      const VectorHandler& XCurr, 
					      const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering SphericalHingeJoint::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indici del vincolo */
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);   
   }
   
   F = Vec3(XCurr, iFirstReactionIndex+1);
   
   const Vec3& x1(pNode1->GetXCurr());
   const Vec3& x2(pNode2->GetXCurr());
   
   Vec3 dTmp1(pNode1->GetRCurr()*d1);
   Vec3 dTmp2(pNode2->GetRCurr()*d2);
   
   WorkVec.Sub(1, F);
   WorkVec.Sub(4, dTmp1.Cross(F));
   WorkVec.Add(7, F);
   WorkVec.Add(10, dTmp2.Cross(F));
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   ASSERT(dCoef != 0.);
   WorkVec.Add(13, (x1+dTmp1-x2-dTmp2)/dCoef);

   if (fc) {
		// Fold = F;
		// Qold = Q;
		bool ChangeJac(false);
		doublereal modF = F.Norm();
		doublereal modFold = Fold.Norm();
		reset_Q = false;
		compute_Q = true;
		// std::cout << "Force: " << F << std::endl;
		if (modFold < preF && modF > preF) {
			// std::cout << "Condition 1" << std::endl;
			// reset_Q = true;
			// std::cout << "reset_Q: " << reset_Q << "; compute_Q: " << compute_Q << std::endl;
		}
		if (modFold > preF && modF > preF && F.Dot(Fold)/modF/modFold < 5.E-1) {
			// std::cout << "Condition 2 " << F.Dot(Fold)/modF/modFold << std::endl;
			// reset_Q = true;
			// std::cout << "reset_Q: " << reset_Q << "; compute_Q: " << compute_Q << std::endl;
		}
		// if (reset_Q) {
		// 	const_cast<VectorHandler&>(XCurr).PutCoef(iFirstReactionIndex+NumSelfDof+1, 0.);
		// 	const_cast<VectorHandler&>(XCurr).PutCoef(iFirstReactionIndex+NumSelfDof+2, 0.);
		// 	const_cast<VectorHandler&>(XPrimeCurr).PutCoef(iFirstReactionIndex+NumSelfDof+1, 0.);
		// 	const_cast<VectorHandler&>(XPrimeCurr).PutCoef(iFirstReactionIndex+NumSelfDof+2, 0.);
		// }
		if (modF <= preF / 2.) {
			// std::cout << "Condition 3" << std::endl;
			compute_Q = false;
			// std::cout << "reset_Q: " << reset_Q << "; compute_Q: " << compute_Q << std::endl;
		} else {
			// std::cout << "Call Spherical" << std::endl;
			// std::cout << "F: " << F << std::endl;
			// std::cout << "reset_Q: " << reset_Q << "; compute_Q: " << compute_Q << std::endl;
			SphericalQR(F, Q, !(reset_Q), Qold);
			// std::cout << "Qold: " << Qold << std::endl;
			// std::cout << "Q   : " << Q << std::endl;
		}
		// std::cout << modF << " " << preF << std::endl;
		if (reset_Q) {
			Qold = Q;
		}
		// std::cout << Q << std::endl;
		// std::cout << "----------------------------" << std::endl;
		modF = std::max(modF, preF);
		const Vec3& Omega1(pNode1->GetWCurr());
		const Vec3& Omega2(pNode2->GetWCurr());
		Vec3 Omegar(Omega1 - Omega2);
		d2D v;
		v.x[0] = (-Q.GetCol(1)).Cross(Omegar).Dot(Q.GetCol(2))*r;
		v.x[1] = (-Q.GetCol(1)).Cross(Omegar).Dot(Q.GetCol(3))*r;
		try {
			fc->AssRes(WorkVec,12+NumSelfDof,iFirstReactionIndex+NumSelfDof,modF,v,XCurr,XPrimeCurr);
		}
		catch (Elem::ChangedEquationStructure& err) {
			ChangeJac = true;
		}
		d2D f = fc->fc();
		d2D shc = Sh_c->Sh_c(f, modF, v);
		// std::cout << "WorkVec  before f" << std::endl;
		// std::cout << WorkVec << std::endl;
		// if (compute_Q) {
		// FIXME: F.Norm() or modF?
			Ffrict1 = -Q.GetCol(2)*modF*f.x[0];
			Ffrict2 = -Q.GetCol(3)*modF*f.x[1];
			// std::cout << "..............." << std::endl;
			// std::cout << "Omega1: " << Omega1 << std::endl;
			// std::cout << "Omega2: " << Omega2 << std::endl;
			// std::cout << "Omegar: " << Omegar << std::endl;
			// std::cout << "F: " << F << std::endl;
			// std::cout << "modF: " << modF << std::endl;
			// std::cout << "t0: " << Q.GetCol(1) << std::endl;
			// std::cout << "t1: " << Q.GetCol(2) << std::endl;
			// std::cout << "t2: " << Q.GetCol(3) << std::endl;
			// std::cout << "v.x[0]: " << v.x[0] << std::endl;
			// std::cout << "v.x[1]: " << v.x[1] << std::endl;
			// std::cout << "Ffrict1: " << Ffrict1 << "; f1: " << fc->fc().x[0] << std::endl;
			// std::cout << "Ffrict2: " << Ffrict2 << "; f2: " << fc->fc().x[1] << std::endl;
			// std::cout << "Ftot: " << F + Ffrict1 + Ffrict2 << std::endl;
			WorkVec.Sub(1, Ffrict1);
			WorkVec.Sub(1, Ffrict2);
			WorkVec.Add(7, Ffrict1);
			WorkVec.Add(7, Ffrict2);
			WorkVec.Sub(4, dTmp1.Cross(Ffrict1)); /* Sfrutto  F/\d = -d/\F */
			WorkVec.Sub(4, dTmp1.Cross(Ffrict2)); /* Sfrutto  F/\d = -d/\F */
			// std::cerr << "F: " << F << std::endl;
			// std::cerr << "Ffrict1: " << Ffrict1 << std::endl;
			// std::cerr << "Ffrict2: " << Ffrict2 << std::endl;
			// std::cerr << "modF: " << modF << std::endl;
			// std::cerr << "Q1: " << Q.GetCol(1) << std::endl;
			// std::cerr << "Q2: " << Q.GetCol(2) << std::endl;
			// std::cerr << "Q3: " << Q.GetCol(3) << std::endl;
			// std::cerr << "m1 0: " << dTmp2.Cross(Ffrict1) << std::endl;
			// std::cerr << "m2 0: " << dTmp2.Cross(Ffrict2) << std::endl;
			WorkVec.Add(10, dTmp2.Cross(Ffrict1));
			WorkVec.Add(10, dTmp2.Cross(Ffrict2));
		// }
		// else {
		// 	Ffrict1 = Zero3;
		// 	Ffrict2 = Zero3;
		// }
		M1 = Q.GetCol(1).Cross(Q.GetCol(2))* shc.x[0] * r * modF;
		M2 = Q.GetCol(1).Cross(Q.GetCol(3))* shc.x[1] * r * modF;
		// std::cout << "M1: " << M1 << std::endl;
		// std::cout << "M2: " << M2 << std::endl;
		// 	std::cout << "Omegar: " << Omegar << std::endl;
		// 	std::cout << "v: " << v.x[0] << " " << v.x[1] << std::endl;
		// 	std::cout << "F: " << F << std::endl;
		// 	std::cout << "Ffrict1: " << Ffrict1 << "; f1: " << f.x[0] << std::endl;
		// 	std::cout << "Ffrict2: " << Ffrict2 << "; f2: " << f.x[1] << std::endl;
		// 	std::cout << "Ftot: " << F + Ffrict1 + Ffrict2 << std::endl;
		// 	std::cout << "Q: " << Q << std::endl;
		// 	std::cout << "Qold: " << Qold << std::endl;
		// std::cout << "Q1: " << Q.GetCol(1) << std::endl;
		// std::cout << "Q2: " << Q.GetCol(2) << std::endl;
		// std::cout << "Q3: " << Q.GetCol(3) << std::endl;
		// std::cout << "M1: " << M1 << std::endl;
		// std::cout << "M2: " << M2 << std::endl;
		WorkVec.Sub(4, M1);
		WorkVec.Sub(4, M2);
		// std::cerr << "M1: " << M1 << std::endl;
		// std::cerr << "M2: " << M2 << std::endl;
		WorkVec.Add(10, M1);
		WorkVec.Add(10, M2);

		if (ChangeJac) {
			throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
		}
	}
	// std::cout << WorkVec << std::endl;
   return WorkVec;
}

			    
DofOrder::Order SphericalHingeJoint::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i >=0 and i < iGetNumDof(), 
		"INDEX ERROR in SphericalHingeJoint::GetEqType");
	if (i<NumSelfDof) {
		return DofOrder::ALGEBRAIC;
	} else {
		return fc->GetEqType(i-NumSelfDof);
	}
}

void
SphericalHingeJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			OutputPrepare_int("Spherical hinge", OH);

			Var_Phi = OH.CreateRotationVar(m_sOutputNameBase, "", od, 
				"relative orientation, in joint reference frame");
			if (fc) {
				Var_MFR = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "MFR",
						OutputHandler::Dimensions::Moment,
						"Overall fricton moment");

				Var_n = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "n",
						OutputHandler::Dimensions::Length,
						"direction n ");
				Var_t1 = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "t1",
						OutputHandler::Dimensions::Length,
						"direction t1 ");
				Var_t2 = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "t2",
						OutputHandler::Dimensions::Length,
						"direction t2 ");

				Var_fc1 = OH.CreateVar<doublereal>(m_sOutputNameBase + "." "fc1",
						OutputHandler::Dimensions::Dimensionless,
						"friction model specific data: friction coefficient in direction t1");
				Var_fc2 = OH.CreateVar<doublereal>(m_sOutputNameBase + "." "fc2",
						OutputHandler::Dimensions::Dimensionless,
						"friction model specific data: friction coefficient in direction t2");

				Var_Fn = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "Fn",
						OutputHandler::Dimensions::Force,
						"Frictionless reaction force in direction n");
				Var_F1 = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "F1",
						OutputHandler::Dimensions::Force,
						"Friction force in direction t1");
				Var_F2 = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "F2",
						OutputHandler::Dimensions::Force,
						"Friction force in direction t2");


			}
		}
#endif // USE_NETCDF
	}
}

/* Output (da mettere a punto) */
void
SphericalHingeJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		Mat3x3 R1Tmp(pNode1->GetRCurr()*R1h);
		Mat3x3 RTmp(R1Tmp.MulTM(pNode2->GetRCurr()*R2h));
		Vec3 E;
		switch (od) {
		case EULER_123:
			E = MatR2EulerAngles123(RTmp)*dRaDegr;
			break;

		case EULER_313:
			E = MatR2EulerAngles313(RTmp)*dRaDegr;
			break;

		case EULER_321:
			E = MatR2EulerAngles321(RTmp)*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			E = RotManip::VecRot(RTmp);
			break;

		case ORIENTATION_MATRIX:
			break;

		default:
			/* impossible */
			break;
		}

		Vec3 Ftot = F+Ffrict1+Ffrict2;
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, (R1Tmp.MulTV(Ftot)), Zero3, Ftot, Zero3);
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				OH.WriteNcVar(Var_Phi, E);
				break;

			case ORIENTATION_MATRIX:
				OH.WriteNcVar(Var_Phi, RTmp);
				break;

			default:
				/* impossible */
				break;
			}
			if (fc) {
				Vec3 Mtot = M1 + M2;
				OH.WriteNcVar(Var_MFR, Mtot);
				OH.WriteNcVar(Var_n,  Q.GetCol(1));
				OH.WriteNcVar(Var_t1, Q.GetCol(2));
				OH.WriteNcVar(Var_t2, Q.GetCol(3));
				OH.WriteNcVar(Var_fc1, fc->fc().x[0]);
				OH.WriteNcVar(Var_fc2, fc->fc().x[1]);
				OH.WriteNcVar(Var_Fn, F);
				OH.WriteNcVar(Var_F1, Ffrict1);
				OH.WriteNcVar(Var_F2, Ffrict2);
			}
		}
#endif // USE_NETCDF

		if (OH.UseText(OutputHandler::JOINTS)) {
			std::ostream &of = Joint::Output(OH.Joints(), "SphericalHinge", GetLabel(),
				R1Tmp.MulTV(Ftot), Zero3, Ftot, Zero3)
				<< " ";

			switch (od) {
				case EULER_123:
				case EULER_313:
				case EULER_321:
				case ORIENTATION_VECTOR:
					OH.Joints() << E;
					break;

				case ORIENTATION_MATRIX:
					OH.Joints() << RTmp;
					break;

				default:
					/* impossible */
					break;
			}
			if(fc) {
				of << " " << M1 + M2;
				of << " " <<  -Q.GetCol(1);
				of << " " <<  Q.GetCol(2);
				of << " " <<  Q.GetCol(3);
				of << " " <<  fc->fc().x[0];
				of << " " <<  fc->fc().x[1];
				of << " " <<  F;
				of << " " <<  Ffrict1;
				of << " " << Ffrict2;
			}
			of << std::endl;
		}
	}
}

void
SphericalHingeJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh == 0) {
				continue;
			}

			if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
				const Mat3x3& R1(pNode1->GetRCurr());
				Vec3 dTmp2(pNode2->GetRCurr()*d2);
   
				d1 = R1.MulTV(pNode2->GetXCurr() + dTmp2 - pNode1->GetXCurr());

			} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
				const Mat3x3& R2(pNode2->GetRCurr());
				Vec3 dTmp1(pNode1->GetRCurr()*d1);
   
				d2 = R2.MulTV(pNode1->GetXCurr() + dTmp1 - pNode2->GetXCurr());

			} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
				/* TODO */
			}
		}
	}
	if (fc) {
		fc->SetValue(pDM, X, XP, ph, iGetFirstIndex() + NumSelfDof);
	}
}

Hint *
SphericalHingeJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /*}*/, STRLENOF("offset{" /*}*/)) == 0) {
		s += STRLENOF("offset{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '2':
			return new Joint::OffsetHint<2>;
		}
	} else if (fc) {
		return fc->ParseHint(pDM, s);
	}

	return 0;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
SphericalHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				   const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering SphericalHingeJoint::InitialAssJac()" << std::endl);

   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Equazioni: vedi joints.dvi */
    
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;

   /* Setto gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
      WM.PutColIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
      WM.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   /* Matrici identita' */
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      /* Contributo di forza all'equazione della forza, nodo 1 */
      WM.PutCoef(iCnt, 24+iCnt, 1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 1 */
      WM.PutCoef(6+iCnt, 27+iCnt, 1.);
      
      /* Contributo di forza all'equazione della forza, nodo 2 */
      WM.PutCoef(12+iCnt, 24+iCnt, -1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 2 */
      WM.PutCoef(18+iCnt, 27+iCnt, -1.);
      
      /* Equazione di vincolo, nodo 1 */
      WM.PutCoef(24+iCnt, iCnt, -1.);
      
      /* Derivata dell'equazione di vincolo, nodo 1 */
      WM.PutCoef(27+iCnt, 6+iCnt, -1.);
      
      /* Equazione di vincolo, nodo 2 */
      WM.PutCoef(24+iCnt, 12+iCnt, 1.);
      
      /* Derivata dell'equazione di vincolo, nodo 2 */
      WM.PutCoef(27+iCnt, 18+iCnt, 1.);
   }
   
   /* Recupera i dati */
   const Mat3x3& R1(pNode1->GetRRef());
   const Mat3x3& R2(pNode2->GetRRef());
   const Vec3& Omega1(pNode1->GetWRef());
   const Vec3& Omega2(pNode2->GetWRef());
   /* F e' stata aggiornata da InitialAssRes */
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   
   /* Distanza nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);

   /* Matrici F/\d1/\, -F/\d2/\ */
   Mat3x3 FWedged1Wedge(MatCrossCross, F, d1Tmp);
   Mat3x3 FWedged2Wedge(MatCrossCross, F, -d2Tmp);
   
   /* Matrici (omega1/\d1)/\, -(omega2/\d2)/\ */
   Mat3x3 O1Wedged1Wedge(MatCross, Omega1.Cross(d1Tmp));
   Mat3x3 O2Wedged2Wedge(MatCross, d2Tmp.Cross(Omega2));
   
   /* Equazione di momento, nodo 1 */
   WM.Add(4, 4, FWedged1Wedge);
   WM.Add(4, 25, Mat3x3(MatCross, d1Tmp));
   
   /* Equazione di momento, nodo 2 */
   WM.Add(16, 16, FWedged2Wedge);
   WM.Sub(16, 25, Mat3x3(MatCross, d2Tmp));
   
   /* Derivata dell'equazione di momento, nodo 1 */
   WM.Add(10, 4, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega1))*Mat3x3(MatCross, d1Tmp));
   WM.Add(10, 10, FWedged1Wedge);
   WM.Add(10, 25, O1Wedged1Wedge);
   WM.Add(10, 28, Mat3x3(MatCross, d1Tmp));
   
   /* Derivata dell'equazione di momento, nodo 2 */
   WM.Sub(22, 16, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega2))*Mat3x3(MatCross, d2Tmp));
   WM.Add(22, 22, FWedged2Wedge);
   WM.Add(22, 25, O2Wedged2Wedge);
   WM.Sub(22, 28, Mat3x3(MatCross, d2Tmp));
      
   /* Equazione di vincolo */
   WM.Add(25, 4, Mat3x3(MatCross, d1Tmp));
   WM.Sub(25, 16, Mat3x3(MatCross, d2Tmp));
   
   /* Derivata dell'equazione di vincolo */
   WM.Add(28, 4, O1Wedged1Wedge);
   WM.Add(28, 10, Mat3x3(MatCross, d1Tmp));
   WM.Add(28, 16, O2Wedged2Wedge);
   WM.Sub(28, 22, Mat3x3(MatCross, d2Tmp));
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
SphericalHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering SphericalHingeJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
   
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
      WorkVec.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati */
   const Vec3& x1(pNode1->GetXCurr());
   const Vec3& x2(pNode2->GetXCurr());
   const Vec3& v1(pNode1->GetVCurr());
   const Vec3& v2(pNode2->GetVCurr());
   const Mat3x3& R1(pNode1->GetRCurr());
   const Mat3x3& R2(pNode2->GetRCurr());
   const Vec3& Omega1(pNode1->GetWCurr());
   const Vec3& Omega2(pNode2->GetWCurr());
   F = Vec3(XCurr, iFirstReactionIndex+1);
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   
   /* Distanza nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);

   /* Vettori omega1/\d1, -omega2/\d2 */
   Vec3 O1Wedged1(Omega1.Cross(d1Tmp));
   Vec3 O2Wedged2(Omega2.Cross(d2Tmp));
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(d1Tmp)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Derivate delle equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(7, FPrime);
   WorkVec.Add(10, FPrime.Cross(d1Tmp)-O1Wedged1.Cross(F));
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(13, F);
   WorkVec.Add(16, d2Tmp.Cross(F)); 
   
   /* Derivate delle equazioni di equilibrio, nodo 2 */
   WorkVec.Add(19, FPrime);
   WorkVec.Add(22, d2Tmp.Cross(FPrime)+O2Wedged2.Cross(F));
   
   /* Equazione di vincolo */
   WorkVec.Add(25, x1+d1Tmp-x2-d2Tmp);
   
   /* Deivata dell'equazione di vincolo */
   WorkVec.Add(28, v1+O1Wedged1-v2-O2Wedged2);
      
   return WorkVec;
}

const OutputHandler::Dimensions
SphericalHingeJoint::GetEquationDimension(integer index) const {
	// DOF == 3
	OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::Length;
			break;
		case 2:
			dimension = OutputHandler::Dimensions::Length;
			break;
      case 3:
			dimension = OutputHandler::Dimensions::Length;
			break;
	  default:
			if (fc) {
				index -= NumSelfDof;
				integer iFCDofs = fc->iGetNumDof();
				if (iFCDofs > 0) {
					/* TODO */
					/* not sure this is handled correctly */
					dimension = fc->GetEquationDimension(index);
				}
			} else {
				dimension = OutputHandler::Dimensions::UnknownDimension;
			}
			break;
	}

	return dimension;
}

/* SphericalHingeJoint - end */


/* PinJoint - begin */
/* Costruttore non banale */
PinJoint::PinJoint(unsigned int uL, const DofOwner* pDO,	       
		   const StructNode* pN,
		   const Vec3& X0Tmp, const Vec3& dTmp,
		   const OrientationDescription& od_a,
		   flag fOut)
: Joint(uL, pDO, fOut), pNode(pN), X0(X0Tmp), d(dTmp), 
F(Zero3),
od(od_a)
{
   NO_OP;
}


/* Distruttore banale */
PinJoint::~PinJoint(void)
{
   NO_OP;
};


/* Contributo al file di restart */
std::ostream& PinJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", pin, "
     << pNode->GetLabel() << ", reference, node, ",
     d.Write(out, ", ") << ", reference, global, ",
     X0.Write(out, ", ") << ';' << std::endl;
   
   return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
PinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		 doublereal dCoef,
		 const VectorHandler& /* XCurr */ ,
		 const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PinJoint::AssJac()" << std::endl);
      
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(27, 0);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   const Mat3x3& R(pNode->GetRRef());
   Vec3 dTmp(R*d);
      
   /* 
    * L'equazione di vincolo afferma che il punto in cui si trova la
    * cerniera deve essere fissato:
    *      x + d = x0
    * 
    * con: d = R * d_0
    * 
    * La forza e' data dalla reazione vincolare F, nel sistema globale
    * La coppia dovuta all'eccentricita' e' data rispettivamente da:
    *     d /\ F
    *
    * 
    *       x      g         F
    * Q1 |  0      0         I   | | x |   | -F          |
    * G1 |  0      cF/\d1/\  d/\ | | g |   | -d/\F       |
    * F  |  I      d/\       0   | | F |   |  (x+d-x0)/c |
    * 
    * con d = R*d_0, c = dCoef
    */

   /* termini di reazione sul nodo */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutItem(iCnt, iFirstMomentumIndex+iCnt, 
		  iFirstReactionIndex+iCnt, 1.);
   }   
   WM.PutCross(4, iFirstMomentumIndex+3,
		iFirstReactionIndex, dTmp);
      
   /* Nota: F, la reazione vincolare, e' stata aggiornata da AssRes */
   
   /* Termini diagonali del tipo: c*F/\d/\Delta_g 
    * nota: la forza e' gia' moltiplicata per dCoef */      
   WM.PutMat3x3(10, iFirstMomentumIndex+3,
		 iFirstPositionIndex+3, Mat3x3(MatCrossCross, F*dCoef, dTmp));

   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* termini di vincolo dovuti al nodo 1 */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutItem(18+iCnt, iFirstReactionIndex+iCnt, 
		  iFirstPositionIndex+iCnt, -1.);
   }
   WM.PutCross(22, iFirstReactionIndex,
		iFirstPositionIndex+3, dTmp);
         
   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& PinJoint::AssRes(SubVectorHandler& WorkVec,
					      doublereal dCoef,
					      const VectorHandler& XCurr,
					      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PinJoint::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
      
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }
     
   
   /* Indici del vincolo */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);
   }

   F = Vec3(XCurr, iFirstReactionIndex+1);
   
   const Vec3& x(pNode->GetXCurr());
   const Mat3x3& R(pNode->GetRCurr());
   
   Vec3 dTmp(R*d);
   
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(dTmp)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   ASSERT(dCoef != 0.);
   WorkVec.Add(7, (x+dTmp-X0)/dCoef);

   return WorkVec;
}

DofOrder::Order PinJoint::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i >=0 and i < iGetNumDof(), 
		"INDEX ERROR in PinJoint::GetEqType");
	return DofOrder::ALGEBRAIC;
}


void
PinJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			OutputPrepare_int("Spherical pin", OH);

			Var_Phi = OH.CreateRotationVar(m_sOutputNameBase, "", od, 
				"relative orientation, in joint reference frame");
		}
#endif // USE_NETCDF
	}
}


/* Output (da mettere a punto) */
void PinJoint::Output(OutputHandler& OH) const
{
#if 0
	if (bToBeOutput()) {
		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "Pin", GetLabel(), F, Zero3, F, Zero3) 
				<< " " << MatR2EulerAngles(pNode->GetRCurr())*dRaDegr << std::endl;
			OH.Joints() << std::endl;
		}
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, F, Zero3, F, Zero3);
			OH.WriteNcVar(Var_Phi, MatR2EulerAngles(pNode->GetRCurr())*dRaDegr);
		}
#endif // USE_NETCDF
	}
#endif 

	if (bToBeOutput()) {
		Vec3 E;
		switch (od) {
		case EULER_123:
			E = MatR2EulerAngles123(pNode->GetRCurr())*dRaDegr;
			break;

		case EULER_313:
			E = MatR2EulerAngles313(pNode->GetRCurr())*dRaDegr;
			break;

		case EULER_321:
			E = MatR2EulerAngles321(pNode->GetRCurr())*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			E = RotManip::VecRot(pNode->GetRCurr());
			break;

		case ORIENTATION_MATRIX:
			break;

		default:
			/* impossible */
			break;
		}
      
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, pNode->GetRCurr().MulTV(F), Zero3, F, Zero3);
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				OH.WriteNcVar(Var_Phi, E);
				break;

			case ORIENTATION_MATRIX:
				OH.WriteNcVar(Var_Phi, pNode->GetRCurr());
				break;

			default:
				/* impossible */
				break;
			}
		}
#endif // USE_NETCDF

		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "SphericalHinge", GetLabel(),
				pNode->GetRCurr().MulTV(F), Zero3, F, Zero3)
				<< " ";

			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				OH.Joints() << E;
				break;

			case ORIENTATION_MATRIX:
				OH.Joints() << pNode->GetRCurr();
				break;

			default:
				/* impossible */
				break;
			}

			OH.Joints() << std::endl;
		}
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
PinJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				   const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering PinJoint::InitialAssJac()" << std::endl);

   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Equazioni: vedi joints.dvi */
    
   /* Indici */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;

   /* Setto gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
      WM.PutColIndex(6+iCnt, iFirstVelocityIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   /* Matrici identita' */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      /* Contributo di forza all'equazione della forza */
      WM.PutCoef(iCnt, 12+iCnt, 1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza */
      WM.PutCoef(6+iCnt, 15+iCnt, 1.);
      
      /* Equazione di vincolo */
      WM.PutCoef(12+iCnt, iCnt, -1.);
      
      /* Derivata dell'equazione di vincolo */
      WM.PutCoef(15+iCnt, 6+iCnt, -1.);
   }
   
   /* Recupera i dati */
   const Mat3x3& R(pNode->GetRRef());
   const Vec3& Omega(pNode->GetWRef());
   /* F e' stata aggiornata da InitialAssRes */
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   
   /* Distanza nel sistema globale */
   Vec3 dTmp(R*d);

   /* Matrici F/\d/\ */
   Mat3x3 FWedgedWedge(MatCrossCross, F, dTmp);
   
   /* Matrici (omega/\d)/\ */
   Mat3x3 OWedgedWedge(MatCross, Omega.Cross(dTmp));
   
   /* Equazione di momento */
   WM.Add(4, 4, FWedgedWedge);
   WM.Add(4, 13, Mat3x3(MatCross, dTmp));
   
   /* Derivata dell'equazione di momento */
   WM.Add(10, 4, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega))*Mat3x3(MatCross, dTmp));
   WM.Add(10, 10, FWedgedWedge);
   WM.Add(10, 13, OWedgedWedge);
   WM.Add(10, 16, Mat3x3(MatCross, dTmp));
   
   /* Equazione di vincolo */
   WM.Add(13, 4, Mat3x3(MatCross, dTmp));
   
   /* Derivata dell'equazione di vincolo */
   WM.Add(16, 4, OWedgedWedge);
   WM.Add(16, 10, Mat3x3(MatCross, dTmp));
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
PinJoint::InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering PinJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
   
   /* Indici */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati */
   const Vec3& x(pNode->GetXCurr());
   const Vec3& v(pNode->GetVCurr());
   const Mat3x3& R(pNode->GetRCurr());
   const Vec3& Omega(pNode->GetWCurr());
   F = Vec3(XCurr, iFirstReactionIndex+1);
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   
   /* Distanza nel sistema globale */
   Vec3 dTmp(R*d);

   /* Vettori omega/\d */
   Vec3 OWedged(Omega.Cross(dTmp));
   
   /* Equazioni di equilibrio */
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(dTmp)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Derivate delle equazioni di equilibrio */
   WorkVec.Sub(7, FPrime);
   WorkVec.Add(10, FPrime.Cross(dTmp)-OWedged.Cross(F));
   
   /* Equazione di vincolo */
   WorkVec.Add(13, x+dTmp-X0);
   
   /* Derivata dell'equazione di vincolo */
   WorkVec.Add(16, v+OWedged);
      
   return WorkVec;
}

const OutputHandler::Dimensions
PinJoint::GetEquationDimension(integer index) const {
	// DOF == 3
   OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::Length;
			break;
		case 2:
			dimension = OutputHandler::Dimensions::Length;
			break;
      case 3:
			dimension = OutputHandler::Dimensions::Length;
			break;
	}

	return dimension;
}

std::ostream&
PinJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": " <<
			"position constraints" << std::endl;

	return out;
}
/* PinJoint - end */
