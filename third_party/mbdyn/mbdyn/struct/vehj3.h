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

#ifndef VEHJ3_H
#define VEHJ3_H

#include "joint.h"
#include "constltp.h"


/* DeformableJoint - begin */

class DeformableJoint :
public Joint {
protected:
	const StructNode* pNode1;
	const StructNode* pNode2;
	mutable Vec3 tilde_f1;
	mutable Vec3 tilde_f2;
	mutable Mat3x3 tilde_R1h;
	mutable Mat3x3 tilde_R2h;

	OrientationDescription od;

	// tilde_d, tilde_ThetaCurr
	Vec6 tilde_k;

	// tilde_dPrime, tilde_Omega
	Vec6 tilde_kPrime;

#ifdef USE_NETCDF
	MBDynNcVar Var_tilde_d;
	MBDynNcVar Var_tilde_dPrime;
	MBDynNcVar Var_d;
	MBDynNcVar Var_dPrime;
	MBDynNcVar Var_Phi;
	MBDynNcVar Var_Omega;
#endif // USE_NETCDF

	bool bFirstRes;

	ConstitutiveLaw6D* pDC;

	Vec3 d1, d2;
	Vec3 d1Prime, d2Prime;
	Mat3x3 R1h;

	Vec6 F;

	// uses F, d1, d2
	void
	AssMatCommon(FullSubMatrixHandler& WM,
		doublereal dCoef);

	// uses d1, d2
	void
	AssMatElastic(FullSubMatrixHandler& WM,
		doublereal dCoef, const Mat6x6& FDE);

	// uses d1, d2, d1Prime, d2Prime
	void
	AssMatViscous(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef, const Mat6x6& FDEPrime);

	virtual void
	AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef) = 0;

	virtual void
	AssVec(SubVectorHandler& WorkVec) = 0;

public:
	/* Costruttore non banale */
	DeformableJoint(unsigned int uL,
		const DofOwner* pDO,
		ConstitutiveLaw6D*const pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut);

	/* Distruttore */
	virtual ~DeformableJoint(void);

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const override {
		return Joint::DEFORMABLEJOINT;
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

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0) override;

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const override;
	         
	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const override;
 
	/* Tipo di DeformableJoint */
	virtual ConstLawType::Type GetConstLawType(void) const = 0;

	virtual unsigned int iGetNumDof(void) const override{
		return 0;
	};

	virtual DofOrder::Order GetDofType(unsigned int /* i */ ) const override {
		return DofOrder::UNKNOWN;
	};

	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
		*piNumRows = 12;
		*piNumCols = 12;
	};

	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int iGetInitialNumDof(void) const override {
		return 0;
	};

	virtual unsigned int iGetNumPrivData(void) const override;
	virtual unsigned int iGetPrivDataIdx(const char *s) const override;
	virtual doublereal dGetPrivData(unsigned int i) const override;

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const override {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr) override;

	/* assemblaggio jacobiano */
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr) override;

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr) override;

	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS) override;

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr) override;

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const override;
};

/* DeformableJoint - end */


/* ElasticJoint - begin */

class ElasticJoint : public DeformableJoint {
protected:
	Vec3 ThetaRef;

	Mat6x6 FDE;

	void AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef) override;
	void AssVec(SubVectorHandler& WorkVec) override;

public:
	ElasticJoint(unsigned int uL,
		const DofOwner* pDO,
		ConstitutiveLaw6D*const pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut);

	~ElasticJoint(void);

	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const override {
		return ConstLawType::ELASTIC;
	};

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP) override;

	/* assemblaggio jacobiano */
	using Elem::AssMats;
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr) override;

	/* Inverse Dynamics Jacobian matrix assembly */
	using DeformableJoint::AssJac;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) override;

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP) override;

	/* Inverse Dynamics update */
	void Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS) override;

	virtual void AfterConvergence(const VectorHandler& X,
			const VectorHandler& XP,
			const VectorHandler& XPP) override;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
		*piNumRows = 12;
		*piNumCols = 12;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) override;

};

/* ElasticJoint - end */

/* ElasticJointInv - begin */

class ElasticJointInv : public DeformableJoint {
protected:
	Vec3 ThetaRef;

	Mat6x6 FDE;

	void AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef) override;
	void AssVec(SubVectorHandler& WorkVec) override;

public:
	ElasticJointInv(unsigned int uL,
		const DofOwner* pDO,
		ConstitutiveLaw6D*const pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut);

	~ElasticJointInv(void);

	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const override {
		return ConstLawType::ELASTIC;
	};

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP) override;

	/* assemblaggio jacobiano */
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr) override;

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP) override;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
		*piNumRows = 12;
		*piNumCols = 12;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) override;

};

/* ElasticJoint - end */

/* ViscousJoint - begin */

class ViscousJoint : public DeformableJoint {
protected:
	Mat6x6 FDEPrime;

	using Elem::AssMats;
	void AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef) override;
	void AssVec(SubVectorHandler& WorkVec) override;

public:
	ViscousJoint(unsigned int uL,
		const DofOwner* pDO,
		ConstitutiveLaw6D*const pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1,
		const Mat3x3& tilde_R2,
		const OrientationDescription& od,
		flag fOut);

	~ViscousJoint(void);

	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const override {
		return ConstLawType::VISCOUS;
	};

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP) override;

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP) override;

	/* Aggiorna le deformazioni ecc. */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
		*piNumRows = 12;
		*piNumCols = 12;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) override;

};

/* ViscousJoint - end */


/* ViscoElasticJoint - begin */

class ViscoElasticJoint
: public DeformableJoint {
protected:
	Vec3 ThetaRef;

	Mat6x6 FDE;
	Mat6x6 FDEPrime;

	using Elem::AssMats;
	void AssMats(FullSubMatrixHandler& WorkMatA,
		FullSubMatrixHandler& WorkMatB,
		doublereal dCoef) override;
	void AssVec(SubVectorHandler& WorkVec) override;

public:
	ViscoElasticJoint(unsigned int uL,
		const DofOwner* pDO,
		ConstitutiveLaw6D*const pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1,
		const Mat3x3& tilde_R2,
		const OrientationDescription& od,
		flag fOut);

	~ViscoElasticJoint(void);

	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const override {
		return ConstLawType::VISCOELASTIC;
	};

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP) override;

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP) override;

	/* Aggiorna le deformazioni ecc. */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
		*piNumRows = 12;
		*piNumCols = 12;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) override;

};

/* ViscoElasticJoint - end */


#endif /* VEHJ3_H */

