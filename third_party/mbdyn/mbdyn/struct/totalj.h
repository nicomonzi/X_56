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

#ifndef TOTALJ_H
#define TOTALJ_H

#include "joint.h"
#include "constltp.h"

#if 0
// friction will be dealt with later
#include "friction.h"
#endif

extern unsigned int total_iGetPrivDataIdx(const char *s_in);

/* TotalJoint - begin */

class TotalJoint : public Joint {
protected:
        const StructNode* pNode1;
        const StructNode* pNode2;
        Vec3 f1;
        Mat3x3 R1h;
        Mat3x3 R1hr;
        Vec3 f2;
        Mat3x3 R2h;
        Mat3x3 R2hr;
        bool bPosActive[3];
        bool bRotActive[3];

        bool bVelActive[3];
        bool bAgvActive[3];	/* Agv stands for AnGular Velocity */

        TplDriveOwner<Vec3> XDrv;
        TplDriveOwner<Vec3> XPDrv;
        TplDriveOwner<Vec3> XPPDrv;

        TplDriveOwner<Vec3> ThetaDrv;
        TplDriveOwner<Vec3> OmegaDrv;
        TplDriveOwner<Vec3> OmegaPDrv;

protected:
        unsigned int nConstraints;
        unsigned int nPosConstraints;
        unsigned int nRotConstraints;
        unsigned int nVelConstraints;
        unsigned int nAgvConstraints;

        unsigned int iPosIncid[3];
        unsigned int iRotIncid[3];
        unsigned int iVelIncid[3];
        unsigned int iAgvIncid[3];

        unsigned int iPosEqIndex[3];
        unsigned int iRotEqIndex[3];
        unsigned int iVelEqIndex[3];
        unsigned int iAgvEqIndex[3];

        Vec3 tilde_f1;

#ifdef USE_NETCDF
        MBDynNcVar Var_X;
        MBDynNcVar Var_Phi;
        MBDynNcVar Var_V;
        MBDynNcVar Var_Omega;
#endif // USE_NETCDF

protected:
        mutable Vec3 M;
        mutable Vec3 F;
        mutable Vec3 ThetaDelta;
        mutable Vec3 ThetaDeltaPrev;

        Vec3 ThetaDeltaRemnant;
        Vec3 ThetaDeltaTrue;

public:
        /* Constructor */
        TotalJoint(unsigned int uL, const DofOwner *pDO,
                bool bPos[3], bool bVel[3],
                TplDriveCaller<Vec3> *const pDCPos[3],
                bool bRot[3], bool bAgv[3],
                TplDriveCaller<Vec3> *const pDCRot[3],
                const StructNode* pN1,
                const Vec3& f1Tmp, const Mat3x3& R1hTmp, const Mat3x3& R1hrTmp,
                const StructNode* pN2,
                const Vec3& f2Tmp, const Mat3x3& R2hTmp, const Mat3x3& R2hrTmp,
                flag fOut);

        /* Destructor */
        ~TotalJoint(void);

        /* Contributo al file di restart */
        virtual std::ostream& Restart(std::ostream& out) const override;
        virtual void Restart(RestartData& oData, RestartData::RestartAction eAction) override;

        /* Tipo di Joint */
        virtual Joint::Type
        GetJointType(void) const override {
                return Joint::TOTALJOINT;
        }

        virtual unsigned int
        iGetNumDof(void) const override {
                return nConstraints;
        }

        virtual std::ostream&
        DescribeDof(std::ostream& out,
                const char *prefix = "",
                bool bInitial = false) const override;

        virtual void
        DescribeDof(std::vector<std::string>& desc,
                bool bInitial = false,
                int i = -1) const override;

        virtual std::ostream&
        DescribeEq(std::ostream& out,
                const char *prefix = "",
                bool bInitial = false) const override;

        virtual void
        DescribeEq(std::vector<std::string>& desc,
                bool bInitial = false,
                int i = -1) const override;

        DofOrder::Order
        GetDofType(unsigned int i) const override {
                ASSERT(i >= 0 && i < nConstraints);
                return DofOrder::ALGEBRAIC;
        }

        virtual void
        SetValue(DataManager *pDM,
                VectorHandler& X, VectorHandler& XP,
                SimulationEntity::Hints *ph = 0) override;

        virtual Hint *
        ParseHint(DataManager *pDM, const char *s) const override;

        virtual void
        AfterConvergence(const VectorHandler& X,
                const VectorHandler& XP) override;

        /* Inverse Dynamics: */
        virtual void
        AfterConvergence(const VectorHandler& X,
                const VectorHandler& XP, const VectorHandler& XPP) override;

        void
        WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
                *piNumRows = 12 + nConstraints ;
                *piNumCols = *piNumRows;
        }

        VariableSubMatrixHandler&
        AssJac(VariableSubMatrixHandler& WorkMat,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr) override;

        SubVectorHandler&
        AssRes(SubVectorHandler& WorkVec,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr) override;

        /* inverse dynamics capable element */
        virtual bool bInverseDynamics(void) const override;

        /* Inverse Dynamics Jacobian matrix assembly */
        VariableSubMatrixHandler&
        AssJac(VariableSubMatrixHandler& WorkMat,
                const VectorHandler& XCurr) override;

        /* Inverse Dynamics residual assembly */
        SubVectorHandler&
        AssRes(SubVectorHandler& WorkVec,
                const VectorHandler& XCurr,
                const VectorHandler&  XPrimeCurr,
                const VectorHandler&  XPrimePrimeCurr,
                InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS) override;

        /* Inverse Dynamics update */
        void Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS) override;

        DofOrder::Order GetEqType(unsigned int i) const override;

        void OutputPrepare(OutputHandler &OH) override;
        void Output(OutputHandler& OH) const override;

        /* funzioni usate nell'assemblaggio iniziale */

        virtual unsigned int
        iGetInitialNumDof(void) const override {
                return  2*nConstraints;
        }

        virtual void
        InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
                *piNumCols = *piNumRows = 24 + 2*nConstraints;
        }

        /* Contributo allo jacobiano durante l'assemblaggio iniziale */
        virtual VariableSubMatrixHandler&
        InitialAssJac(VariableSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr) override;

        /* Contributo al residuo durante l'assemblaggio iniziale */
        virtual SubVectorHandler&
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr) override;

        /* Dati privati */
        virtual unsigned int iGetNumPrivData(void) const override;
        virtual unsigned int iGetPrivDataIdx(const char *s) const override;
        virtual doublereal dGetPrivData(unsigned int i) const override;

        /* *******PER IL SOLUTORE PARALLELO******** */
        /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
         * utile per l'assemblaggio della matrice di connessione fra i dofs */
        virtual void
        GetConnectedNodes(std::vector<const Node *>& connectedNodes) const override {
                connectedNodes.resize(2);
                connectedNodes[0] = pNode1;
                connectedNodes[1] = pNode2;
        }
        /* ************************************************ */

        /* returns the dimension of the component */
        const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const override;
};

/* TotalJoint - end */

/* TotalPinJoint - begin */

class TotalPinJoint : public Joint {
protected:
        const StructNode* pNode;
        Vec3 Xc;
        Mat3x3 Rch;
        Mat3x3 Rchr;
        Vec3 tilde_fn;
        Mat3x3 tilde_Rnh;
        Mat3x3 tilde_Rnhr;
        bool bPosActive[3];
        bool bRotActive[3];
        bool bVelActive[3];
        bool bAgvActive[3];	/* Agv stands for AnGular Velocity */

        Mat3x3 RchT;
        Vec3 tilde_Xc;
        Mat3x3 RchrT;

        TplDriveOwner<Vec3> XDrv;
        TplDriveOwner<Vec3> XPDrv;
        TplDriveOwner<Vec3> XPPDrv;

        TplDriveOwner<Vec3> ThetaDrv;
        TplDriveOwner<Vec3> OmegaDrv;
        TplDriveOwner<Vec3> OmegaPDrv;

        unsigned int nConstraints;
        unsigned int nPosConstraints;
        unsigned int nRotConstraints;
        unsigned int nVelConstraints;
        unsigned int nAgvConstraints;

        unsigned int iPosIncid[3];
        unsigned int iRotIncid[3];
        unsigned int iVelIncid[3];
        unsigned int iAgvIncid[3];

        unsigned int iPosEqIndex[3];
        unsigned int iRotEqIndex[3];
        unsigned int iVelEqIndex[3];
        unsigned int iAgvEqIndex[3];

#ifdef USE_NETCDF
        MBDynNcVar Var_X;
        MBDynNcVar Var_Phi;
        MBDynNcVar Var_V;
        MBDynNcVar Var_Omega;
#endif // USE_NETCDF

protected:
        mutable Vec3 M;
        mutable Vec3 F;
        mutable Vec3 ThetaDelta;
        mutable Vec3 ThetaDeltaPrev;

        Vec3 ThetaDeltaRemnant;
        Vec3 ThetaDeltaTrue;

public:
        /* Constructor */
        TotalPinJoint(unsigned int uL, const DofOwner *pDO,
                bool bPos[3], bool bVel[3],
                TplDriveCaller<Vec3> *const pDCPos[3],
                bool bRot[3], bool bAgv[3],
                TplDriveCaller<Vec3> *const pDCRot[3],
                const Vec3& XcTmp, const Mat3x3& RchTmp, const Mat3x3& RchrTmp,
                const StructNode* pN,
                const Vec3& fnTmp, const Mat3x3& RnhTmp, const Mat3x3& RnhrTmp,
                flag fOut);

        /* Destructor */
        ~TotalPinJoint(void);

        /* Contributo al file di restart */
        virtual std::ostream& Restart(std::ostream& out) const override;

        virtual void Restart(RestartData& oData, RestartData::RestartAction eAction) override;

        /* Tipo di Joint */
        virtual Joint::Type
        GetJointType(void) const override {
                return Joint::TOTALPINJOINT;
        }

        virtual unsigned int
        iGetNumDof(void) const override {
                return nConstraints;
        }

        virtual std::ostream&
        DescribeDof(std::ostream& out,
                const char *prefix = "",
                bool bInitial = false) const override;

        virtual void
        DescribeDof(std::vector<std::string>& desc,
                bool bInitial = false,
                int i = -1) const override;

        virtual std::ostream&
        DescribeEq(std::ostream& out,
                const char *prefix = "",
                bool bInitial = false) const override;

        virtual void
        DescribeEq(std::vector<std::string>& desc,
                bool bInitial = false,
                int i = -1) const override;

        DofOrder::Order
        GetDofType(unsigned int i) const override {
                ASSERT(i >= 0 && i < nConstraints);
                return DofOrder::ALGEBRAIC;
        }

        virtual void
        SetValue(DataManager *pDM,
                VectorHandler& X, VectorHandler& XP,
                SimulationEntity::Hints *ph = 0) override;

        virtual Hint *
        ParseHint(DataManager *pDM, const char *s) const override;

        virtual void
        AfterConvergence(const VectorHandler& X,
                         const VectorHandler& XP) override;

        virtual void
        AfterConvergence(const VectorHandler& X,
                         const VectorHandler& XP, const VectorHandler& XPP) override;

        void
        WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
                *piNumRows = 6 + nConstraints ;
                *piNumCols = *piNumRows;
        }

        VariableSubMatrixHandler&
        AssJac(VariableSubMatrixHandler& WorkMat,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr) override;

        SubVectorHandler&
        AssRes(SubVectorHandler& WorkVec,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr) override;

        /* inverse dynamics capable element */
        virtual bool bInverseDynamics(void) const override;

        /* inverse dynamics Jacobian matrix assembly */
        VariableSubMatrixHandler&
        AssJac(VariableSubMatrixHandler& WorkMat,
               const VectorHandler& XCurr) override;

        /* inverse dynamics residual assembly */
        SubVectorHandler&
        AssRes(SubVectorHandler& WorkVec,
                const VectorHandler& XCurr,
                const VectorHandler&  XPrimeCurr,
                const VectorHandler&  XPrimePrimeCurr,
                InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS) override;

        /* Inverse Dynamics update */
        virtual void Update(const VectorHandler& XCurr,
                            InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS) override;


        DofOrder::Order GetEqType(unsigned int i) const override;

        void OutputPrepare(OutputHandler &OH) override;
        void Output(OutputHandler& OH) const override;

        /* funzioni usate nell'assemblaggio iniziale */

        virtual unsigned int
        iGetInitialNumDof(void) const override {
                return  2*nConstraints;
        }
        virtual void
        InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
                *piNumCols = *piNumRows = 12 + 2*nConstraints;
        }

        /* Contributo allo jacobiano durante l'assemblaggio iniziale */
        virtual VariableSubMatrixHandler&
        InitialAssJac(VariableSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr) override;

        /* Contributo al residuo durante l'assemblaggio iniziale */
        virtual SubVectorHandler&
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr) override;

        /* Dati privati */
        virtual unsigned int iGetNumPrivData(void) const override;
        virtual unsigned int iGetPrivDataIdx(const char *s) const override;
        virtual doublereal dGetPrivData(unsigned int i) const override;

        /* *******PER IL SOLUTORE PARALLELO******** */
        /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
         * utile per l'assemblaggio della matrice di connessione fra i dofs */
        virtual void
        GetConnectedNodes(std::vector<const Node *>& connectedNodes) const override {
                connectedNodes.resize(1);
                connectedNodes[0] = pNode;
        }
        /* ************************************************ */

        /* returns the dimension of the component */
        const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const override;
};

/* TotalPinJoint - end */




/****************************************/
/****************************************/
/*              TOTAL FORCE	        */
/****************************************/
/****************************************/

#include "force.h"

/* Total Force: begin */
class TotalForce : public Force	{
private:
        const StructNode* pNode1;
        const StructNode* pNode2;
        Vec3 f1;
        Mat3x3 R1h;
        Mat3x3 R1hr;
        Vec3 f2;
        Mat3x3 R2h;
        Mat3x3 R2hr;

        TplDriveOwner<Vec3> FDrv;

        TplDriveOwner<Vec3> MDrv;

        mutable Vec3 M;
        mutable Vec3 F;

public:
        TotalForce(unsigned int uL,
                TplDriveCaller<Vec3> *const pDCForce,
                TplDriveCaller<Vec3> *const pDCCouple,
                const StructNode* pN1,
                const Vec3& f1Tmp, const Mat3x3& R1hTmp, const Mat3x3& R1hrTmp,
                const StructNode* pN2,
                const Vec3& f2Tmp, const Mat3x3& R2hTmp, const Mat3x3& R2hrTmp,
                flag fOut);

        ~TotalForce(void) {
                NO_OP;
        }

        /* Force Type */
        virtual Force::Type GetForceType(void) const override {
                return Force::TOTALINTERNALFORCE;
        }

        virtual std::ostream& Restart(std::ostream& out) const override;

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
                *piNumRows = 12;
                *piNumCols = 12;
        }

        virtual VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr) override;

        virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr) override;

        /* inverse dynamics capable element */
        virtual bool bInverseDynamics(void) const override;

        /* Inverse Dynamics*/
        virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
                        const VectorHandler& /* XCurr */ ,
                        const VectorHandler& /* XPrimeCurr */ ,
                        const VectorHandler& /* XPrimePrimeCurr */ ,
                        InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS) override;

        virtual void Output(OutputHandler& OH) const override;

        virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
                *piNumRows = 24;
                *piNumCols = 24;
        }

        virtual VariableSubMatrixHandler&
        InitialAssJac(VariableSubMatrixHandler& WorkMat,
                const VectorHandler& XCurr) override;

        virtual SubVectorHandler&
        InitialAssRes(SubVectorHandler& WorkVec,
                const VectorHandler& XCurr) override;
};

/* Total Force: end */

#endif // TOTALJ_H
