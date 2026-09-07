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

/* Deformable hinges */


#ifndef VEHJ_H
#define VEHJ_H

#include "joint.h"
#include "constltp.h"

extern const char* psConstLawNames[];

/* DeformableHingeJoint - begin */

class DeformableHingeJoint :
public Joint {
protected:
	const StructNode* pNode1;
	const StructNode* pNode2;
	mutable Mat3x3 tilde_R1h;
	mutable Mat3x3 tilde_R2h;

	OrientationDescription od;

	ConstitutiveLaw3D* pDC;
private:
#ifdef USE_NETCDF
	MBDynNcVar Var_Phi;
	MBDynNcVar Var_Omega;
#endif // USE_NETCDF

protected:
	bool bFirstRes;

	Vec3 M;

	Mat3x3 MDE;
	Mat3x3 MDEPrime;

	/* for invariant stuff */
	Mat3x3	hat_I;
	Mat3x3	hat_IT;

	/* Jacobian matrix helpers */
	virtual void
	AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef);

	void
	AssMatMInv(FullSubMatrixHandler& WMA, doublereal dCoef);

	void
	AssMatMDE(FullSubMatrixHandler& WMA, doublereal dCoef);

	virtual void
	AssMatMDEPrime(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB, doublereal dCoef);

	void
	AssMatMDEPrimeInv(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB, doublereal dCoef);

	/* output helper */
	void OutputInv(OutputHandler& OH) const;

	/* priv data helper */
	doublereal
	dGetPrivDataInv(unsigned int i) const;
	virtual void AfterPredictHelper(void) = 0;

public:
	/* Costruttore non banale */
	DeformableHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		ConstitutiveLaw3D*const pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut);

	/* Distruttore */
	virtual ~DeformableHingeJoint(void);

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const override {
		return Joint::DEFORMABLEHINGE;
	};
    
	/* Deformable element */
	virtual bool bIsDeformable() const override {
		return true;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const override;

        virtual void Restart(RestartData& oData, RestartData::RestartAction eAction) override;

	void OutputPrepare(OutputHandler& OH) override;
	virtual void Output(OutputHandler& OH) const override;

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP) override;

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0) override;
	virtual void SetInitialValue(VectorHandler& /* X */ ) override;

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const override;
	         
	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const = 0;

	virtual unsigned int iGetNumDof(void) const override {
		return 0;
	};

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
		*piNumRows = 6;
		*piNumCols = 6;
	};

	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const override;
 
	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int iGetInitialNumDof(void) const override {
		return 0;
	};

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const override {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */

	virtual unsigned int iGetNumPrivData(void) const override;
	virtual unsigned int iGetPrivDataIdx(const char *s) const override;
	virtual doublereal dGetPrivData(unsigned int i) const override;

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const override;
};

/* DeformableHingeJoint - end */


/* ElasticHingeJoint - begin */

class ElasticHingeJoint : public DeformableHingeJoint {
protected:
	Vec3 ThetaRef;
	Vec3 ThetaCurr;

	virtual void AfterPredictHelper(void);
	virtual void AssMat(FullSubMatrixHandler& WM, doublereal dCoef);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ElasticHingeJoint(unsigned int uL,
			const DofOwner* pDO,
			ConstitutiveLaw3D*const pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	virtual ~ElasticHingeJoint(void);

	virtual void
	AfterPredict(VectorHandler& X, VectorHandler& XP);

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/* assemblaggio jacobiano */
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
			VariableSubMatrixHandler& WorkMatB,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics Jacobian matrix assembly */
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr,
		const VectorHandler& XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	/* Inverse Dynamics update */
	void Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	virtual void AfterConvergence(const VectorHandler& X,
			const VectorHandler& XP,
			const VectorHandler& XPP);

	virtual void InitialWorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const {
		*piNumRows = 6;
		*piNumCols = 6;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

#if 0
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableHingeJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};

	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableHingeJoint::dGetPrivData(i);
	};
#endif
};

/* ElasticHingeJoint - end */


/* ElasticHingeJointInv - begin */

class ElasticHingeJointInv : public ElasticHingeJoint {
protected:
	virtual void
	AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef);

	/* AssMatMDE is OK as MDE is updated fine by AfterPredict();
	 * AssMatMDEPrime is not needed */

	virtual void AfterPredictHelper(void);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ElasticHingeJointInv(unsigned int uL,
			const DofOwner* pDO,
			ConstitutiveLaw3D*const pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	virtual ~ElasticHingeJointInv(void);

	virtual void Output(OutputHandler& OH) const;

	virtual doublereal dGetPrivData(unsigned int i) const;

};

/* ElasticHingeJointInv - end */


/* ViscousHingeJoint - begin */

class ViscousHingeJoint : public DeformableHingeJoint {
protected:
	Vec3 Omega;

	virtual void AfterPredictHelper(void);
	virtual void AssMats(FullSubMatrixHandler& WMA,
			FullSubMatrixHandler& WMB,
			doublereal dCoef);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscousHingeJoint(unsigned int uL,
			const DofOwner* pDO,
			ConstitutiveLaw3D*const pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	virtual ~ViscousHingeJoint(void);

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOUS;
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/* assemblaggio jacobiano */
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
			VariableSubMatrixHandler& WorkMatB,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	virtual void InitialWorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const  {
		*piNumRows = 6;
		*piNumCols = 12;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

#if 0
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableHingeJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};

	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableHingeJoint::dGetPrivData(i);
	};
#endif
};

/* ViscousHingeJoint - end */


/* ViscousHingeJointInv - begin */

class ViscousHingeJointInv : public ViscousHingeJoint {
protected:
	/* AssMatMDEPrime is not needed */
	virtual void
	AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef);
	virtual void
	AssMatMDEPrime(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB, doublereal dCoef);

	virtual void AfterPredictHelper(void);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscousHingeJointInv(unsigned int uL,
			const DofOwner* pDO,
			ConstitutiveLaw3D*const pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	virtual ~ViscousHingeJointInv(void);

	virtual void Output(OutputHandler& OH) const;

	virtual doublereal dGetPrivData(unsigned int i) const;
};

/* ViscousHingeJointInv - end */


/* ViscoElasticHingeJoint - begin */

class ViscoElasticHingeJoint
: public DeformableHingeJoint {
protected:
	Vec3 ThetaRef;
	Vec3 ThetaCurr;

	Vec3 Omega;

	virtual void AssMats(FullSubMatrixHandler& WMA,
			FullSubMatrixHandler& WMB,
			doublereal dCoef);

	virtual void AfterPredictHelper(void);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscoElasticHingeJoint(unsigned int uL,
			const DofOwner* pDO,
			ConstitutiveLaw3D*const pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	~ViscoElasticHingeJoint(void);

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
			VariableSubMatrixHandler& WorkMatB,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	virtual void InitialWorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const {
		*piNumRows = 6;
		*piNumCols = 12;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

#if 0
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableHingeJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};

	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableHingeJoint::dGetPrivData(i);
	};
#endif
};

/* ViscoElasticHingeJoint - end */

/* ViscoElasticHingeJointInv - begin */

class ViscoElasticHingeJointInv
: public ViscoElasticHingeJoint {
protected:
	/* AssMatMDEPrime is not needed */
	virtual void
	AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef);
	virtual void
	AssMatMDEPrime(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB, doublereal dCoef);

	virtual void AfterPredictHelper(void);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscoElasticHingeJointInv(unsigned int uL,
			const DofOwner* pDO,
			ConstitutiveLaw3D*const pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	~ViscoElasticHingeJointInv(void);

	virtual void Output(OutputHandler& OH) const;

	virtual doublereal dGetPrivData(unsigned int i) const;
};

/* ViscoElasticHingeJointInv - end */


/* InvAngularCLR - begin */

struct InvAngularCLR : public ConstitutiveLawRead<Vec3, Mat3x3> {
	virtual ConstitutiveLaw<Vec3, Mat3x3> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);
};

/* InvAngularCLR - end */

#endif /* VEHJ_H */

