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

/* Driven elements:
 * elements that are used depending on the (boolean) value
 * of a driver. Example: a driven joint is assembled only
 * if the driver is true, otherwise there is no joint and
 * the reaction unknowns are set to zero
 */

#ifndef DRIVEN_H
#define DRIVEN_H

#include "nestedelem.h"
#include "drive.h"

#include "except.h"

class DrivenElem : 
	public NestedElem, protected DriveOwner {
protected:
	DataManager *pDM;
	SimulationEntity::Hints *pHints;
	bool bActive;

#ifdef USE_NETCDF
	MBDynNcVar Var_status;
#endif // USE_NETCDF

public:
	DrivenElem(DataManager *pDM, const DriveCaller* pDC, bool b_active,
			const Elem* pE, SimulationEntity::Hints *ph = 0);
	~DrivenElem(void);

	virtual bool bIsActive(void) const;

	virtual void OutputPrepare(OutputHandler& OH) override;
	virtual void Output(OutputHandler& OH) const override;

	virtual void SetValue(DataManager *pdm,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0) override;

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const override;

        virtual void Restart(RestartData& oData, RestartData::RestartAction eAction) override;

	/* funzioni proprie */

	/*
	 * Elaborazione vettori e dati prima e dopo la predizione
	 * per MultiStepIntegrator */
	virtual void BeforePredict(VectorHandler& X,
		VectorHandler& XP,
		std::deque<VectorHandler*>& qXPr,
		std::deque<VectorHandler*>& qXPPr) const override;

	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP) override;

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr) override;

	virtual void AfterConvergence(const VectorHandler& X,
     			const VectorHandler& XP) override;

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
	    		const VectorHandler& XCurr,
	    		const VectorHandler& XPrimeCurr) override;

     	virtual void AssMats(VariableSubMatrixHandler& WorkMatA,
 			VariableSubMatrixHandler& WorkMatB,
 			const VectorHandler& XCurr,
 			const VectorHandler& XPrimeCurr) override;

	/* assemblaggio residuo */
     	virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr) override;

	/*
	 * Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData()
	 */
	virtual doublereal dGetPrivData(unsigned int i) const override;

	/* Inverse Dynamics: */
	virtual void Update(const VectorHandler& XCurr,
			InverseDynamics::Order iOrder) override;

	/* inverse dynamics Jacobian matrix assembly */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) override;

	/* inverse dynamics residual assembly */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr,
		const VectorHandler& XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS) override;

	/* Inverse Dynamics: */
	virtual void AfterConvergence(const VectorHandler& X,
     		const VectorHandler& XP, const VectorHandler& XPP) override;

	/* InitialAssemblyElem */
public:
	virtual unsigned int iGetInitialNumDof(void) const override;

	/* Dimensione del workspace durante l'assemblaggio iniziale. Occorre tener
	 * conto del numero di dof che l'elemento definisce in questa fase e dei
	 * dof dei nodi che vengono utilizzati. Sono considerati dof indipendenti
	 * la posizione e la velocita' dei nodi */
	virtual void InitialWorkSpaceDim(integer* piNumRows,
		integer* piNumCols) const override;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) override;

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr) override;

	/* ElemGravityOwner */
protected:
	virtual Vec3 GetS_int(void) const override;
	virtual Mat3x3 GetJ_int(void) const override;

	virtual Vec3 GetB_int(void) const override;

	// NOTE: gravity owners must provide the momenta moment
	// with respect to the origin of the global reference frame!
	virtual Vec3 GetG_int(void) const override;

public:
	virtual doublereal dGetM(void) const override;
	Vec3 GetS(void) const;
	Mat3x3 GetJ(void) const;

	/* ElemDofOwner */
public:
	virtual void SetInitialValue(VectorHandler& X) override;
};

#endif /* DRIVEN_H */

