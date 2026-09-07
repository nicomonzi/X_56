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
  Copyright (C) 2024(-2024) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <cstring>
#include <ctime>

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include <userelem.h>

#include "module-mpc.h"
#include "strnodead.h"
#include "sp_matvecass.h"

using namespace sp_grad;


class MPC: public UserDefinedElem
{
public:
        using Elem::AssRes;
        using Elem::AssJac;
        using UserDefinedElem::InitialAssRes;
        using UserDefinedElem::InitialAssJac;

        MPC(unsigned uLabel, const DofOwner *pDO,
                       DataManager* pDM, MBDynParser& HP);
        virtual ~MPC(void);
        virtual unsigned int iGetNumDof(void) const override;
        virtual DofOrder::Order GetDofType(unsigned int i) const override;
        virtual DofOrder::Order GetEqType(unsigned int i) const override;
        virtual std::ostream& DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const override;
        virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const override;
        virtual unsigned int iGetNumPrivData(void) const override;
        virtual unsigned int iGetPrivDataIdx(const char *s) const override;
        virtual doublereal dGetPrivData(unsigned int i) const override;
        virtual void Output(OutputHandler& OH) const override;
        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;
        virtual VariableSubMatrixHandler&
        AssJac(VariableSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr) override;
        virtual void
        AssJac(VectorHandler& JacY,
               const VectorHandler& Y,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               VariableSubMatrixHandler& WorkMat) override;
        SubVectorHandler&
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr) override;
        template <typename T>
        inline void
        AssRes(SpGradientAssVec<T>& WorkVec,
               doublereal dCoef,
               const SpGradientVectorHandler<T>& XCurr,
               const SpGradientVectorHandler<T>& XPrimeCurr,
               enum SpFunctionCall func);
        int iGetNumConnectedNodes(void) const;
        void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const override;
        void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
                      SimulationEntity::Hints *ph) override;
        std::ostream& Restart(std::ostream& out) const override;
        virtual unsigned int iGetInitialNumDof(void) const override;
        virtual void
        InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;
        VariableSubMatrixHandler&
        InitialAssJac(VariableSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr) override;
        SubVectorHandler&
        InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr) override;
        template <typename T>
        inline void
        InitialAssRes(SpGradientAssVec<T>& WorkVec,
                      const SpGradientVectorHandler<T>& XCurr,
                      enum SpFunctionCall func);

private:
     void SaveLambda(const sp_grad::SpColVector<doublereal>& lambda) {
          for (size_t i = 0; i < rgConstrData.size(); ++i) {
               rgConstrData[i].dLambda = lambda(i + 1);
          }
     }

     void SaveLambda(const sp_grad::SpColVector<sp_grad::SpGradient>&) {
     }

     void SaveLambda(const sp_grad::SpColVector<sp_grad::GpGradProd>&) {
     }

     template <typename T>
     inline void
     UnivAssRes(SpGradientAssVec<T>& WorkVec,
                doublereal dCoef,
                const SpGradientVectorHandler<T>& XCurr,
                enum SpFunctionCall func);

     struct DofData {
          StructDispNodeAd* pNode;
          Vec3 X0;
     };

     struct ConstraintData {
          doublereal dLambda;
          DriveOwner constraintValue;
     };

     std::vector<DofData> rgDofData;
     std::vector<ConstraintData> rgConstrData;
     sp_grad::SpMatrix<doublereal> C;
};

MPC::MPC(unsigned uLabel_a, const DofOwner *pDO,
        DataManager* pDM, MBDynParser& HP)
:       UserDefinedElem(uLabel_a, pDO)
{
        // help
        if (HP.IsKeyWord("help")) {
                silent_cout(
                        "\n"
                        "Module:        mpc\n"
                        "\n"
                        "	This element implements a MPC constraint for nodal displacements\n"
                        "\n"
                        "       mcp,\n"
                        "       nodes, (integer) <number_of_nodes>,\n"
                        "          <node_1>, ..., <node_N>,\n"
                        "       constraints, (integer) <number_of_constraints>\n"
                        "          (DriveCaller) <constr_1>,\n"
                        "          ...\n"
                        "          (DriveCaller) <constr_N>,\n"
                        "       constraint matrix\n"
                        "          (real) <C11>, (real) <C12>, ...\n"
                        "          (real) <C21>, (real) <C22>, ...\n"
                        "          ...\n"
                        "\n"
                        "\n");

                if (!HP.IsArg()) {
                        /*
                         * Exit quietly if nothing else is provided
                         */
                        throw NoErr(MBDYN_EXCEPT_ARGS);
                }
        }

        if ( !HP.IsKeyWord("nodes") ) {
             silent_cerr("mpc(" << GetLabel() << "): keyword \"nodes\" expected at line " << HP.GetLineData() << "\n");
             throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const integer iNumNodes = HP.GetInt();

        if (iNumNodes < 1) {
             silent_cerr("mpc(" << GetLabel() << "): minimum one node is required at line " << HP.GetLineData() << "\n");
             throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        rgDofData.resize(iNumNodes);

        for (integer i = 1; i <= iNumNodes; ++i) {
             rgDofData[i - 1].pNode = pDM->ReadNode<StructDispNodeAd, Node::STRUCTURAL>(HP);
             rgDofData[i - 1].X0 = rgDofData[i - 1].pNode->GetXCurr();
        }

        if ( !HP.IsKeyWord("constraints") ) {
             silent_cerr("mpc(" << GetLabel() << "): keyword \"constraints\" expected at line " << HP.GetLineData() << "\n");
             throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const integer iNumConstr = HP.GetInt();

        if (iNumConstr < 1) {
             silent_cerr("mpc(" << GetLabel() << "): minimum one constraint is required at line " << HP.GetLineData() << "\n");
             throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        rgConstrData.resize(iNumConstr);

        for (integer i = 1; i <= iNumConstr; ++i) {
             rgConstrData[i - 1].constraintValue.Set(HP.GetDriveCaller());
             rgConstrData[i - 1].dLambda = 0.;
        }

        if ( !HP.IsKeyWord("constraint" "matrix") ) {
             silent_cerr("mpc(" << GetLabel() << "): keyword \"constraint matrix\" expected at line " << HP.GetLineData() << "\n");
             throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        C.ResizeReset(iNumConstr, iNumNodes * 3, 0.);

        for (integer i = 1; i <= iNumConstr; ++i) {
             for (integer j = 1; j <= iNumNodes * 3; ++j) {
                  C(i, j) = HP.GetReal();
             }
        }

        SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

MPC::~MPC(void)
{
        // destroy private data
}

unsigned int MPC::iGetNumDof(void) const
{
     return rgConstrData.size();
}

DofOrder::Order MPC::GetDofType(unsigned int i) const
{
     return DofOrder::ALGEBRAIC;
}

DofOrder::Order MPC::GetEqType(unsigned int i) const
{
     return DofOrder::ALGEBRAIC;
}

std::ostream& MPC::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
        const integer iFirstIndex = iGetFirstIndex();

        out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + rgConstrData.size() << ": reaction forces [lambda[1], lambda[2], ...]\n";

        return out;
}

std::ostream& MPC::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
        const integer iFirstIndex = iGetFirstIndex();

        out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + rgConstrData.size() << ": constraint equations [c[1], c[2], ...]\n";

        return out;
}

unsigned int MPC::iGetNumPrivData(void) const
{
     return rgConstrData.size();
}

unsigned int MPC::iGetPrivDataIdx(const char *s) const
{

     unsigned int idx = 0;

     sscanf(s, "lambda[%u]", &idx);

     if (idx > rgConstrData.size()) {
          idx = 0;
     }

     return idx;
}

doublereal MPC::dGetPrivData(unsigned int i) const
{
     if (i <= 0 || i > rgConstrData.size()) {
          silent_cerr("mpc(" << GetLabel() << "): invalid private data index " << i << "\n");
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     return rgConstrData[i - 1].dLambda;
}

void
MPC::Output(OutputHandler& OH) const
{
        if ( bToBeOutput() )
        {
                if ( OH.UseText(OutputHandler::LOADABLE) )
                {
                     std::ostream& os = OH.Loadable();

                     os << std::setw(8) << GetLabel();

                     for (unsigned i = 0; i < rgConstrData.size(); ++i) {
                          os << ' ' << rgConstrData[i].dLambda;
                     }

                     os << '\n';
                }
        }
}

void
MPC::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = rgConstrData.size() + 3 * rgDofData.size();
     *piNumCols = 0;
}

VariableSubMatrixHandler&
MPC::AssJac(VariableSubMatrixHandler& WorkMat,
                       doublereal dCoef,
                       const VectorHandler& XCurr,
                       const VectorHandler& XPrimeCurr)
{
     using namespace sp_grad;

        SpGradientAssVec<SpGradient>::AssJac(this,
                                             WorkMat.SetSparseGradient(),
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC);

        return WorkMat;
}

void
MPC::AssJac(VectorHandler& JacY,
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
MPC::AssRes(SubVectorHandler& WorkVec,
                       doublereal dCoef,
                       const VectorHandler& XCurr,
                       const VectorHandler& XPrimeCurr)
{
     using namespace sp_grad;

        SpGradientAssVec<doublereal>::AssRes(this,
                                             WorkVec,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_RES);

        return WorkVec;
}

template <typename T>
inline void
MPC::UnivAssRes(SpGradientAssVec<T>& WorkVec,
                doublereal dCoef,
                const SpGradientVectorHandler<T>& XCurr,
                enum SpFunctionCall func) {

     const integer iFirstIndex = iGetFirstIndex();

     using namespace sp_grad;

     SpColVector<T> U(3 * rgDofData.size(), 0);
     SpColVectorA<T, 3> Xi;

     for (size_t i = 0; i < rgDofData.size(); ++i) {
          rgDofData[i].pNode->GetXCurr(Xi, dCoef, func);

          for (index_type j = 1; j <= 3; ++j) {
               U(3 * i + j) = Xi(j) - rgDofData[i].X0(j);
          }
     }

     SpColVector<T> Cu = C * U;

     for (size_t i = 0; i < rgConstrData.size(); ++i) {
          Cu(i + 1) -= rgConstrData[i].constraintValue.dGet();
          Cu(i + 1) /= dCoef;
     }

     for (index_type i = 1; i <= Cu.iGetNumRows(); ++i) {
          WorkVec.AddItem(iFirstIndex + i, Cu(i));
     }

     SpColVector<T> lambda(rgConstrData.size(), 0);

     XCurr.GetVec(iFirstIndex + 1, lambda, 1.); // Note: for algebraic variables dCoef is always one

     SaveLambda(lambda);

     SpColVector<T> CTlambda = Transpose(C) * lambda;

     for (size_t i = 0; i < rgDofData.size(); ++i) {
          for (index_type j = 1; j <= 3; ++j) {
               WorkVec.AddItem(rgDofData[i].pNode->iGetFirstMomentumIndex() + j, CTlambda(i * 3 + j));
          }
     }
}

template <typename T>
inline void
MPC::AssRes(SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const SpGradientVectorHandler<T>& XCurr,
            const SpGradientVectorHandler<T>& XPrimeCurr,
            enum SpFunctionCall func) {

     UnivAssRes(WorkVec, dCoef, XCurr, func);
}


int
MPC::iGetNumConnectedNodes(void) const
{
     return rgDofData.size();
}

void
MPC::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
     connectedNodes.resize(0);
     connectedNodes.reserve(rgDofData.size());

     for (const auto& dof: rgDofData) {
          connectedNodes.push_back(dof.pNode);
     }
}

void
MPC::SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
              SimulationEntity::Hints *ph)
{
     const integer iFirstIndex = iGetFirstIndex();

     for (size_t i = 1; i <= rgConstrData.size(); ++i) {
          X(iFirstIndex + i) = rgConstrData[i - 1].dLambda;
     }
}

std::ostream&
MPC::Restart(std::ostream& out) const
{
        return out;
}

unsigned int
MPC::iGetInitialNumDof(void) const
{
     return rgConstrData.size();
}

void
MPC::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = rgConstrData.size() + 3 * rgDofData.size();
     *piNumCols = 0;
}

VariableSubMatrixHandler&
MPC::InitialAssJac(
        VariableSubMatrixHandler& WorkMat,
        const VectorHandler& XCurr)
{

        SpGradientAssVec<SpGradient>::InitialAssJac(this,
                                                    WorkMat.SetSparseGradient(),
                                                    XCurr,
                                                    INITIAL_ASS_JAC);

        return WorkMat;
}

SubVectorHandler&
MPC::InitialAssRes(
        SubVectorHandler& WorkVec,
        const VectorHandler& XCurr)
{
        SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                    WorkVec,
                                                    XCurr,
                                                    INITIAL_ASS_RES);

        return WorkVec;
}

template <typename T>
inline void
MPC::InitialAssRes(SpGradientAssVec<T>& WorkVec,
                   const SpGradientVectorHandler<T>& XCurr,
                   enum SpFunctionCall func) {
     UnivAssRes(WorkVec, 1., XCurr, func);
}

bool mpc_set(void)
{
        UserDefinedElemRead *rf = new UDERead<MPC>;

        if (!SetUDE("mpc", rf))
        {
                delete rf;
                return false;
        }

        return true;
}

#ifndef STATIC_MODULES

extern "C"
{

int module_init(const char *module_name, void *pdm, void *php)
{
        if (!mpc_set())
        {
                silent_cerr("mpc: "
                            "module_init(" << module_name << ") "
                            "failed" << std::endl);

                return -1;
        }

        return 0;
}

}

#endif // ! STATIC_MODULE
