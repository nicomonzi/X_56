/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2023
 *
 * Marco Morandini	<morandini@aero.polimi.it>
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

/* Copyright (C) 2003 Marco Morandini*/

#ifndef FRICTION2D_H
#define FRICTION2D_H

#include "ScalarFunctions.h"
#include "simentity.h"
#include "elem.h"
#include "JacSubMatrix.h"

struct d2D {
	doublereal x[2];
	doublereal operator*(const d2D y) const {return x[0]*y.x[0] + x[1]*y.x[1];}
	d2D operator-(const d2D y) const {return d2D({x[0] - y.x[0], x[1] - y.x[1]});}
	d2D operator+(const d2D y) const {return d2D({x[0] + y.x[0], x[1] + y.x[1]});}
};

std::ostream& operator << (std::ostream& out, const d2D& y);

/** Base class for friction models
 */
class BasicFriction2D : public SimulationEntity{
public:
	using SimulationEntity::SetValue;
	using SimulationEntity::AfterConvergence;
/*
 * 	unsigned int iGetNumDof(void) const;
 * 	DofOrder::Order GetDofType(unsigned int i) const;
 * 	DofOrder::Order GetEqType (unsigned int i) const;
 * 	void SetValue(VectorHandler&X, VectorHandler&XP);
 * 	void BeforePredict(VectorHandler&,
 * 		VectorHandler&,
 * 		VectorHandler&,
 * 		VectorHandler&) const;
 * 	void AfterPredict(VectorHandler&X, VectorHandler&XP);
 * 	void Update(const VectorHandler&XCurr, const VectorHandler&XPrimeCurr);
 * 	void AfterConvergence(VectorHandler&X, VectorHandler&XP);
 */
/** Set Initial Values
 */
	virtual void SetValue(DataManager *pDM,
			VectorHandler&X, VectorHandler&XP,
			SimulationEntity::Hints *ph = 0,
			const unsigned int solution_startdof = 0);
/** Return last computed friction coefficient
 */
	virtual d2D fc(void) const = 0;
/** Compute self residual and friction coefficient
 */
	virtual void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal F,
		const d2D v,
		const VectorHandler& X,
		const VectorHandler& XP) 
			/*throw(Elem::ChangedEquationStructure)*/ = 0;
/** Compute self jacobian and friction coefficient derivatives
 */
	virtual void AssJac(
		FullSubMatrixHandler& WorkMat,
		ExpandableMatrix& dfc,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal dCoef,
		const doublereal F,
		const d2D v,
		const VectorHandler& X,
		const VectorHandler& XP,
		const ExpandableRowVector& dF,
		const ExpandableMatrix& dv) const = 0;
	virtual void AfterConvergence(
		const doublereal F,
		const d2D v,
		const VectorHandler&X, 
		const VectorHandler&XP,
		const unsigned int solution_startdof) {};
	
	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const = 0;
};

/** Base class for friction shape coefficient
 */
class BasicShapeCoefficient2D {
public:
	virtual ~BasicShapeCoefficient2D(void){};
/** Return last computed shape coefficient
 */
	virtual d2D Sh_c(void) const = 0;
/** Compute the shape coefficient
 */
	virtual d2D Sh_c(
		const d2D f,
		const doublereal F,
		const d2D v) = 0;
/** Compute derivatives of the shape coefficient
 */
	virtual void dSh_c(
		ExpandableMatrix& dShc,
		const d2D f,
		const doublereal F,
		const d2D v,
		const ExpandableMatrix& dfc,
		const ExpandableRowVector& dF,
		const ExpandableMatrix& dv) const = 0;
};

/** A friction model based on 
 * "Dupont Pierre, Hayward, Vincent, Armstrong Brian
 * and Altpeter Friedhelm, Single state elasto-plastic friction models,
 * IEEE Transactions on Automatic Control, scheduled for June 2002"
 */
class ModLugreFriction2D : public BasicFriction2D {
private:
	const doublereal sigma0;
	const doublereal sigma1;
	const doublereal sigma2;
	const doublereal kappa;
	const DifferentiableScalarFunction & fss;
	d2D f;
	doublereal alphatilde(const doublereal zm,
		const doublereal vm) const;
	doublereal alphatilded_zm(const doublereal zm,
		const doublereal vm) const;
	doublereal alphatilded_vm(const doublereal zm,
		const doublereal vm) const;
	doublereal epsilon(const d2D& z, const d2D& v) const;
	doublereal alpha(const d2D& z,
		const d2D& v) const;
	void alphad_v(const d2D& z,
		const d2D& v, ExpandableMatrix& alpha_v) const;
	void alphad_z(const d2D& z,
		const d2D& v, ExpandableRowVector& alpha_z,
		const unsigned int solution_startdof) const;
	const doublereal fs(const doublereal vm) const;
	const doublereal fsd_vm(const doublereal vm) const;
public:
	ModLugreFriction2D(
		const doublereal sigma0,
		const doublereal sigma1,
		const doublereal sigma2,
		const doublereal kappa,
		const BasicScalarFunction *const ff);
	void SetValue(DataManager *pDM,
			VectorHandler&X, VectorHandler&XP,
			SimulationEntity::Hints *ph = 0,
			const unsigned int solution_startdof = 0);
	unsigned int iGetNumDof(void) const;
	virtual std::ostream&
	DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	virtual std::ostream&
	DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	DofOrder::Order GetEqType (unsigned int i) const;
	d2D fc(void) const;
	void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal F,
		const d2D v,
		const VectorHandler& X,
		const VectorHandler& XP) /*throw(Elem::ChangedEquationStructure)*/;
	void AssJac(
		FullSubMatrixHandler& WorkMat,
		ExpandableMatrix& dfc,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal dCoef,
		const doublereal F,
		const d2D v,
		const VectorHandler& X,
		const VectorHandler& XP,
		const ExpandableRowVector& dF,
		const ExpandableMatrix& dv) const;

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

class DiscreteCoulombFriction2D : public BasicFriction2D {
private:
	enum tr_type{
		null,
		from_sticked_to_sliding,
		from_sticking_to_sliding,
		from_sliding_to_sticked,
		from_sliding_to_sticking};
	enum status_type{
		sticked,
		sticking,
		sliding};
	//logical converged_sticked;
	status_type status;
	tr_type transition_type;
	d2D converged_v;
	logical first_iter;
	logical first_switch;
	d2D previous_switch_v;
	d2D current_velocity;
	d2D saved_sliding_velocity;
	d2D saved_sliding_friction;
	doublereal sigma2;
	doublereal vel_ratio;
	d2D current_friction_force;
	doublereal vel_tolerance;

	const DifferentiableScalarFunction & fss;
	d2D f;
	mutable ExpandableMatrix Direction_d;
	bool use_sliding_v;
public:
	DiscreteCoulombFriction2D(
		const BasicScalarFunction *const ff,
		const doublereal s2,
		const doublereal vr,
		const doublereal vt);
	void SetValue(DataManager *pDM,
			VectorHandler&X, VectorHandler&XP,
			SimulationEntity::Hints *ph = 0,
			const unsigned int solution_startdof = 0);
	unsigned int iGetNumDof(void) const;
	virtual std::ostream&
	DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	virtual std::ostream&
	DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	DofOrder::Order GetEqType (unsigned int i) const;
	d2D fc(void) const;
	void AfterConvergence(
		const doublereal F,
		const d2D v,
		const VectorHandler&X,
		const VectorHandler&XP,
		const unsigned int solution_startdof);
	void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal F,
		const d2D v,
		const VectorHandler& X,
		const VectorHandler& XP) /*throw(Elem::ChangedEquationStructure)*/;
	void AssJac(
		FullSubMatrixHandler& WorkMat,
		ExpandableMatrix& dfc,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal dCoef,
		const doublereal F,
		const d2D v,
		const VectorHandler& X,
		const VectorHandler& XP,
		const ExpandableRowVector& dF,
		const ExpandableMatrix& dv) const;

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};



/** Simple shape coefficient: 1.
 */
class SimpleShapeCoefficient2D : public BasicShapeCoefficient2D {
private:
	d2D shc;
public:
	SimpleShapeCoefficient2D(void) : BasicShapeCoefficient2D(), shc({1., 1.}) {};
	virtual d2D Sh_c(void) const;
	virtual d2D Sh_c(
		const d2D f,
		const doublereal F,
		const d2D v);
	virtual void dSh_c(
		ExpandableMatrix& dShc,
		const d2D f,
		const doublereal F,
		const d2D v,
		const ExpandableMatrix& dfc,
		const ExpandableRowVector& dF,
		const ExpandableMatrix& dv) const;
};



//---------------------------------------

BasicFriction2D *const ParseFriction2D(MBDynParser& HP,
	DataManager * pDM);

BasicShapeCoefficient2D *const ParseShapeCoefficient2D(MBDynParser& HP);

#endif /* FRICTION2D_H */

