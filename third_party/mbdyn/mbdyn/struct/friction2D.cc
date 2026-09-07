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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>

#include "mbpar.h"
#include "datamanforward.h"
#include "friction2D.h"
#include "submat.h"

//defined in friction.cc
extern int sign(const doublereal x);

std::ostream& operator << (std::ostream& out, const d2D& y) {
	out << y.x[0] << " " << y.x[1];
	return out;
};

doublereal d2DNorm(const d2D& z) {
	return std::sqrt(z.x[0] * z.x[0] + z.x[1] * z.x[1]);
}

d2D operator*(const doublereal& x, const d2D& y) {
	return d2D({x*y.x[0], x*y.x[1]});
}

doublereal Dot(const d2D& x, const d2D& y) {
	return x.x[0]*y.x[0] + x.x[1]*y.x[1];
}

void d2DNorm_d(const d2D& z, d2D& mdz) {
	doublereal zm = d2DNorm(z);
	if (zm < 1.E-6) {
		mdz.x[0] = sign(z.x[0]);
		mdz.x[1] = sign(z.x[1]);
	} else {
		mdz.x[0] = z.x[0] / zm;
		mdz.x[1] = z.x[1] / zm;
	}
}

d2D d2DDirection(const d2D& z) {
	d2D mdz;
	doublereal zm = d2DNorm(z);
	if (zm < 1.E-6) {
		mdz.x[0] = 0.;
		mdz.x[1] = 0.;
	} else {
		mdz.x[0] = z.x[0] / zm;
		mdz.x[1] = z.x[1] / zm;
	}
	return mdz;
}

void d2DDirection_d(const d2D& z, ExpandableMatrix& mdz) {
	doublereal zm = d2DNorm(z);
	doublereal zm3 = std::pow(zm, 3);
	mdz.ReDim(2, 1);
	mdz.SetBlockDim(1, 2);
	if (zm < 1.E-6) {
		mdz.Set(1., 1, 1, 1);
		mdz.Set(0., 1, 1, 2);
		mdz.Set(0., 2, 1, 1);
		mdz.Set(1., 2, 1, 2);
	} else {
		// mdz.Set(z.x[1]*z.x[1] / zm3, 1, 1, 1);
		// mdz.Set(-z.x[0]*z.x[1] / zm, 1, 1, 2);
		// mdz.Set(-z.x[0]*z.x[1] / zm, 2, 1, 1);
		// mdz.Set(z.x[0]*z.x[0] / zm3, 2, 1, 2);
		mdz.Set(-z.x[0]*z.x[0] / zm3 + 1. / zm, 1, 1, 1);
		mdz.Set(-z.x[0]*z.x[1] / zm3, 1, 1, 2);
		mdz.Set(-z.x[0]*z.x[1] / zm3, 2, 1, 1);
		mdz.Set(-z.x[1]*z.x[1] / zm3 + 1. / zm, 2, 1, 2);
	}
	return;
}

bool operator == (const d2D& x, const d2D& y) {
	return (x.x[0] == y.x[0]) && (x.x[1] == y.x[1]);
}

void
BasicFriction2D::SetValue(DataManager *pDM,
		VectorHandler&X, VectorHandler&XP,
		SimulationEntity::Hints *ph,
		const unsigned int solution_startdof)
{
	NO_OP;
}

ModLugreFriction2D::ModLugreFriction2D(
		const doublereal s0,
		const doublereal s1,
		const doublereal s2,
		const doublereal k,
		const BasicScalarFunction *const ff) : 
sigma0(s0),
sigma1(s1),
sigma2(s2),
kappa(k),
fss(dynamic_cast<const DifferentiableScalarFunction&>(*ff)),
f({0., 0.})
{
	NO_OP;
}

void
ModLugreFriction2D::SetValue(DataManager *pDM,
		VectorHandler&X, 
		VectorHandler&XP,
		SimulationEntity::Hints *ph,
		const unsigned int solution_startdof)
{
	X.PutCoef(solution_startdof+1, f.x[0]/sigma0);
	X.PutCoef(solution_startdof+2, f.x[1]/sigma0);
}

unsigned int ModLugreFriction2D::iGetNumDof(void) const {
	return 2;
};

std::ostream&
ModLugreFriction2D::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out
		<< prefix << "[1]: ModLugreFriction2D state 1" << std::endl
		<< prefix << prefix << "[2]: ModLugreFriction2D state 2" << std::endl;
}

void
ModLugreFriction2D::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(2);
	desc[desc.size()-2] = "ModLugreFriction2D state 1";
	desc[desc.size()-1] = "ModLugreFriction2D state 2";
}

std::ostream&
ModLugreFriction2D::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out
		<< prefix << "[1]: ModLugreFriction2D equation 1" << std::endl
		<< prefix << prefix << "[2]: ModLugreFriction2D equation 2" << std::endl;
}

void
ModLugreFriction2D::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(2);
	desc[0] = "ModLugreFriction2D equation 1";
	desc[1] = "ModLugreFriction2D equation 2";
}

DofOrder::Order ModLugreFriction2D::GetDofType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction2D::GetDofType");
	return DofOrder::DIFFERENTIAL;
};

DofOrder::Order ModLugreFriction2D::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction2D::GetEqType");
	return DofOrder::DIFFERENTIAL;
};

const doublereal ModLugreFriction2D::fs(const doublereal vm) const {
	return fss(vm); // *sign(v); FIXME
};

const doublereal ModLugreFriction2D::fsd_vm(const doublereal vm) const {
	return fss.ComputeDiff(vm); // *sign(v); FIXME
};

doublereal ModLugreFriction2D::alphatilde(const doublereal zm,
	const doublereal vm) const {

	doublereal zss = fs(vm)/sigma0;
	doublereal zba = kappa*zss;

	// std::cout << "zss: " << zss << std::endl;
	// std::cout << "zba: " << zba << std::endl;
	// std::cout << "zm: " << zm << std::endl;
	if (zm <= zba) {
	} else if ((zba <= zm) && (zm <= zss)) {
		// std::cout << "xxx: " << 0.5*std::sin(
		// 			M_PI*sigma0*zm/(fs(vm)*(1.-kappa))-M_PI*(1.+kappa)/(2*(1.-kappa))
		// 		) + 0.5 << std::endl;
		return 0.5*std::sin(
					M_PI*sigma0*zm/(fs(vm)*(1.-kappa))-M_PI*(1.+kappa)/(2*(1.-kappa))
				) + 0.5;
	} else {
		return 1.;
	}
	return 0.;
};

doublereal ModLugreFriction2D::epsilon(const d2D& z,
	const d2D& v) const {
	doublereal zn = d2DNorm(z);
	doublereal vn = d2DNorm(v);
	doublereal eps = 1;
	if (zn > 0. && vn > 0.) {
		eps = (v.x[0]*z.x[0] / zn / vn + v.x[1]*z.x[1] / zn / vn + 1.) / 2.;
	}
	return eps;
};

doublereal ModLugreFriction2D::alpha(const d2D& z,
	const d2D& v) const {
	doublereal zm = d2DNorm(z);
	doublereal vm = d2DNorm(v);
	// doublereal at = alphatilde(zm, vm);
	// std::cout << "alphatilde: " << at << std::endl;
	// std::cout << "epsilon: " << epsilon(z, v) << std::endl;
	return alphatilde(zm, vm) * epsilon(z, v);
}

doublereal ModLugreFriction2D::alphatilded_zm(const doublereal zm,
	const doublereal vm) const {

	doublereal zss = fs(vm)/sigma0;
	doublereal zba = kappa*zss;

	if (zm <= zba) {
		return 0.;
	} else if ((zba <= zm) && (zm <= zss)) {
		return 0.5*M_PI*std::cos(M_PI*sigma0*zm/fs(vm)/(1.-kappa)
			-M_PI*(1.+kappa)/2./(1.-kappa))
			*sigma0/fs(vm)/(1.-kappa);
	} else {
		return 0.;
	}

			// der = 0.5*M_PI*sigma0/fs(v)/(1.-kappa)
			// 	*std::cos(M_PI*sigma0*z/fs(v)/(1.-kappa)
			// 	-M_PI*(1.+kappa)/2./(1.-kappa));


	return 0.;
};

doublereal ModLugreFriction2D::alphatilded_vm(const doublereal zm,
	const doublereal vm) const {

	doublereal zss = fs(vm)/sigma0;
	doublereal zba = kappa*zss;

	if (zm <= zba) {
		return 0.;
	} else if ((zba <= zm) && (zm <= zss)) {
		return -M_PI/2.*std::cos(M_PI*sigma0*zm/fs(vm)/(1.-kappa)
			-M_PI*(1.+kappa)/2./(1.-kappa))
			*sigma0*zm/std::pow(fs(vm),2)/(1.-kappa)
			*fsd_vm(vm);
	} else {
		return 0.;
	}

	// doublereal zss = fs(v)/sigma0;
	// doublereal zba = kappa*zss;
	// doublereal der;
			// der = -M_PI/2*std::cos(M_PI*sigma0*z/fs(v)/(1.-kappa)
			// 	-M_PI*(1.+kappa)/2./(1.-kappa))
			// 	*sigma0*z/std::pow(fs(v),2)/(1.-kappa)
			// 	*fsd(v);


	return 0.;
};

void ModLugreFriction2D::alphad_v(const d2D& z,
	const d2D& v, ExpandableMatrix& alpha_v) const {

	// epsilon = (v.x[0]*z.x[0] + v.x[1]*z.x[1] + 1.) / 2.;
	// alpha = alphatilde * epsilon;

	doublereal zm = d2DNorm(z);
	doublereal vm = d2DNorm(v);
	doublereal eps = epsilon(z, v);
	doublereal alphat = alphatilde(zm, vm);
	doublereal alphatd_vm = alphatilded_vm(zm, vm);
	d2D vmd_v;
	d2DNorm_d(v, vmd_v);
	// d2D der;
	// der.x[0] = alphat * z.x[0];
	// der.x[1] = alphat * z.x[1];
	alpha_v.ReDim(1, 1);
	alpha_v.SetBlockDim(1, 2);
	alpha_v.Set(alphat * z.x[0] / 2. + eps * alphatd_vm * vmd_v.x[0], 1, 1, 1);
	alpha_v.Set(alphat * z.x[1] / 2. + eps * alphatd_vm * vmd_v.x[1], 1, 1, 2);
	return;
};

void ModLugreFriction2D::alphad_z(const d2D& z,
	const d2D& v, ExpandableRowVector& alpha_z, const unsigned int startdof) const {

	doublereal zm = d2DNorm(z);
	doublereal vm = d2DNorm(v);
	doublereal eps = epsilon(z, v);
	doublereal alphat = alphatilde(zm, vm);
	doublereal alphatd_zm = alphatilded_zm(zm, vm);
	d2D zmd_z;
	d2DNorm_d(z, zmd_z);
	// d2D der;
	// der.x[0] = alphat * z.x[0] + eps * alphat_dvm * d2DDirection(z.x[0]);
	// der.x[1] = alphat * z.x[1] + eps * alphat_dvm * d2DDirection(z.x[1]);
	alpha_z.ReDim(2);
	alpha_z.Set(alphat * v.x[0] / 2. + eps * alphatd_zm * zmd_z.x[0], 1, startdof+1);
	alpha_z.Set(alphat * v.x[1] / 2. + eps * alphatd_zm * zmd_z.x[1], 2, startdof+2);

	return;
};

d2D ModLugreFriction2D::fc(void) const {
	return f;
};

void ModLugreFriction2D::AssRes(
	SubVectorHandler& WorkVec,
	const unsigned int startdof,
	const unsigned int solution_startdof,
	const doublereal F,
	const d2D v,
	const VectorHandler& X,
	const VectorHandler& XP) /*throw(Elem::ChangedEquationStructure)*/ {

	d2D z = {X(solution_startdof+1), X(solution_startdof+2)};
	d2D zp = {XP(solution_startdof+1), XP(solution_startdof+2)};

	doublereal vm = d2DNorm(v);
	doublereal fsvm = fs(vm);
	doublereal alph = alpha(z,v);

	f = sigma0*z + sigma1*zp + sigma2*v;
	// f.x[0] = sigma0*z.x[0] + sigma1*zp.x[0] + sigma2*v.x[0];
	// f.x[1] = sigma0*z.x[1] + sigma1*zp.x[1] + sigma2*v.x[1];
	// std::cout << "sigma0: " << sigma0 << std::endl;
	// std::cout << "sigma1: " << sigma1 << std::endl;
	// std::cout << "sigma2: " << sigma2 << std::endl;
	// std::cout << "f: " << f << std::endl;
	// std::cout << "alph: " << alph << std::endl;
	// std::cout << "fsvm: " << fsvm << std::endl;
	// std::cout << "sigma0*z.x[0]: " << sigma0*z.x[0] << std::endl;
	// std::cout << "sigma1*zp.x[0]: " << sigma1*zp.x[0] << std::endl;
	// std::cout << "sigma2*v.x[0]: " << sigma2*v.x[0] << std::endl;
	// std::cout << "z.x[0]:" << z.x[0] << "; zp.x[0]: " << zp.x[0] << "; v.x[0]: " << v.x[0] << std::endl;
	// std::cout << "z.x[1]:" << z.x[1] << "; zp.x[1]: " << zp.x[1] << "; v.x[1]: " << v.x[1] << std::endl;
	// std::cout << "zp: " << zp << std::endl;
	// std::cout << "v: " << v << std::endl;
	// std::cout << "z * alph / fsvm * sigma0: " << alph * z.x[0] / fsvm * sigma0 << " " << alph * z.x[1] / fsvm * sigma0 << std::endl;
	WorkVec.IncCoef(startdof+1, zp.x[0] - v.x[0] + alph * z.x[0] / fsvm * sigma0 * vm);
	WorkVec.IncCoef(startdof+2, zp.x[1] - v.x[1] + alph * z.x[1] / fsvm * sigma0 * vm);
};

void ModLugreFriction2D::AssJac(
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
	const ExpandableMatrix& dv) const {

	d2D z = {X(solution_startdof+1), X(solution_startdof+2)};

	doublereal vm = d2DNorm(v);
	//doublereal zp = XP(solution_startdof+1);

/*
 * 	attrito
 */
 	dfc.ReDim(2, 2);
	dfc.SetBlockDim(1, 2);
	dfc.SetBlockDim(2, 2);
	dfc.SetBlockIdx(1, startdof+1);
	dfc.Set(sigma0*dCoef+sigma1, 1, 1, 1);
	dfc.Set(0., 1, 1, 2);
	dfc.Set(0., 2, 1, 1);
	dfc.Set(sigma0*dCoef+sigma1, 2, 1, 2);
	dfc.Set(sigma2, 1, 2, 1);
	dfc.Set(0., 1, 2, 2);
	dfc.Set(0., 2, 2, 1);
	dfc.Set(sigma2, 2, 2, 2);
	dfc.Link(2, &dv);

/*
 * 	z
 */

	doublereal alph = alpha(z,v);

	ExpandableRowVector alpha_d_z;
	alphad_z(z, v, alpha_d_z, startdof);
	ExpandableMatrix alpha_d_v;
	alphad_v(z, v, alpha_d_v);
	alpha_d_v.Link(1, &dv);

	doublereal fsvm = fs(vm);
	doublereal fsvm2 = fsvm * fsvm;
	doublereal fsvmd_vm = fsd_vm(vm);
	d2D vmd_v;
	d2DNorm_d(v, vmd_v);

	// WorkVec.IncCoef(startdof+1, zp.x[0] - v.x[0] + alph * z.x[0] / fsvm * sigma0);
	// WorkVec.IncCoef(startdof+2, zp.x[1] - v.x[1] + alph * z.x[1] / fsvm * sigma0);

	// -dot(dz) - alpha * sigma0 / fss * dz
	WorkMat.IncCoef(startdof+1, startdof+1, -1. - alph * sigma0 / fsvm * vm * dCoef);
	WorkMat.IncCoef(startdof+2, startdof+2, -1. - alph * sigma0 / fsvm * vm * dCoef);

	ExpandableMatrix deq_dv;
	deq_dv.ReDim(2, 2);

	// -z * sigma0 / fss * alpha_{/v} * dv
	deq_dv.SetBlockDim(1, 1);
	deq_dv.Set(-z.x[0] * sigma0 / fsvm * vm, 1, 1, 1);
	deq_dv.Set(-z.x[1] * sigma0 / fsvm * vm, 2, 1, 1);
	deq_dv.Link(1, &alpha_d_v);

	// dv + z * alpha * sigma0 / fss^2 * fss_{/vm} vm_{/v} * dv
	deq_dv.SetBlockDim(2, 2);
	deq_dv.Set(1. + z.x[0] * alph * sigma0 / fsvm2 * fsvmd_vm * vmd_v.x[0] * vm
					- z.x[0] * alph * sigma0 / fsvm * vmd_v.x[0], 1, 2, 1);
	deq_dv.Set(   + z.x[0] * alph * sigma0 / fsvm2 * fsvmd_vm * vmd_v.x[1] * vm
					- z.x[0] * alph * sigma0 / fsvm * vmd_v.x[1], 1, 2, 2);
	deq_dv.Set(   + z.x[1] * alph * sigma0 / fsvm2 * fsvmd_vm * vmd_v.x[0] * vm
					- z.x[1] * alph * sigma0 / fsvm * vmd_v.x[0], 2, 2, 1);
	deq_dv.Set(1. + z.x[1] * alph * sigma0 / fsvm2 * fsvmd_vm * vmd_v.x[1] * vm
					- z.x[1] * alph * sigma0 / fsvm * vmd_v.x[1], 2, 2, 2);
	deq_dv.Link(2, &dv);

	deq_dv.AddTo(WorkMat, startdof+1);

	// -z sigma0/fss * alpha_{/z} * dz
	ExpandableMatrix deq_dz;
	deq_dz.ReDim(2, 1);
	deq_dz.SetBlockDim(1, 1);
	//deq_dz.SetBlockIdx(1, startdof+1);
	deq_dz.Set(-z.x[0] * sigma0 / fsvm * vm, 1, 1, 1);
	deq_dz.Set(-z.x[1] * sigma0 / fsvm * vm, 2, 1, 1);
	deq_dz.Link(1, &alpha_d_z);
	deq_dz.AddTo(WorkMat, startdof+1, dCoef);


//	std::cout << alphad_z(z,v) << std::endl;
/*
 * 	callback: dfc[] = df/d{F,v,(z+dCoef,zp)}
 */
};

const OutputHandler::Dimensions
ModLugreFriction2D::GetEquationDimension(integer index) const {
	// DOF == 1
	OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::Velocity;
			break;
	}

	return dimension;
}

//-------------------------

DiscreteCoulombFriction2D::DiscreteCoulombFriction2D(
		const BasicScalarFunction *const ff,
		const doublereal s2,
		const doublereal vr,
		const doublereal vt) :
//converged_sticked(true),
status(sticked),
transition_type(null),
converged_v({0., 0.}),
first_iter(true),
first_switch(true),
previous_switch_v({0., 0.}),
current_velocity({0., 0.}),
sigma2(s2),
vel_ratio(vr),
current_friction_force({0., 0.}),
vel_tolerance(vt),
fss(dynamic_cast<const DifferentiableScalarFunction&>(*ff)),
f({0., 0.})
{
	NO_OP;
}

void
DiscreteCoulombFriction2D::SetValue(DataManager *pDM,
		VectorHandler&X,
		VectorHandler&XP,
		SimulationEntity::Hints *ph,
		const unsigned int solution_startdof)
{
	X.PutCoef(solution_startdof+1,f.x[0]);
	X.PutCoef(solution_startdof+2,f.x[1]);
}

unsigned int DiscreteCoulombFriction2D::iGetNumDof(void) const {
	return 2;
};

std::ostream&
DiscreteCoulombFriction2D::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out
		<< prefix << "[1]: DiscreteCoulombFriction2D state 1" << std::endl
		<< prefix << prefix << "[2]: DiscreteCoulombFriction2D state 2" << std::endl;
}

void
DiscreteCoulombFriction2D::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(2);
	desc[desc.size()-2] = "DiscreteCoulombFriction2D state 1";
	desc[desc.size()-1] = "DiscreteCoulombFriction2D state 2";
}

std::ostream&
DiscreteCoulombFriction2D::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out
		<< prefix << "[1]: DiscreteCoulombFriction2D equation 1" << std::endl
		<< prefix << prefix << "[2]: DiscreteCoulombFriction2D equation 2" << std::endl;
}

void
DiscreteCoulombFriction2D::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(2);
	desc[desc.size()-2] = "DiscreteCoulombFriction2D equation 1";
	desc[desc.size()-1] = "DiscreteCoulombFriction2D equation 2";
}

DofOrder::Order DiscreteCoulombFriction2D::GetDofType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in DiscreteCoulombFriction2D::GetDofType");
	return DofOrder::ALGEBRAIC;
};

DofOrder::Order DiscreteCoulombFriction2D::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in DiscreteCoulombFriction2D::GetEqType");
	return DofOrder::DIFFERENTIAL;
};

d2D DiscreteCoulombFriction2D::fc(void) const {
	return current_friction_force;
};

void DiscreteCoulombFriction2D::AfterConvergence(
	const doublereal F,
	const d2D v,
	const VectorHandler&X,
	const VectorHandler&XP,
	const unsigned int solution_startdof) {
	f.x[0] = X(solution_startdof+1);
	f.x[1] = X(solution_startdof+2);
	converged_v = v;
	current_velocity = v;
	previous_switch_v = v;
	transition_type = null;
	first_iter = true;
	first_switch = true;
	if (status == sticking) {
		status = sticked;
	} else if (status == sliding) {
	} else {
	}
};


void DiscreteCoulombFriction2D::AssRes(
	SubVectorHandler& WorkVec,
	const unsigned int startdof,
	const unsigned int solution_startdof,
	const doublereal F,
	const d2D v,
	const VectorHandler& X,
	const VectorHandler& XP)  {
	f.x[0] = X(solution_startdof+1);
	f.x[1] = X(solution_startdof+2);
	transition_type = null;
	use_sliding_v = false;
	// std::cerr << "f: " << f << std::endl;
	// std::cerr << "fss(0): " << fss(0) << std::endl;
	// std::cerr << "d2DNorm(f): " << d2DNorm(f) << std::endl;
	if (d2DNorm(f)-fss(0) > 1.0E-6*fss(0)) {
		//unconditionally switch to sliding
		if (status == sticked) {
			transition_type = from_sticked_to_sliding;
		} else if (status == sticking) {
			transition_type = from_sticking_to_sliding;
		} else if (status == sliding) {
			//do nothing
		} else {
			silent_cerr("DiscreteCoulombFriction2D::AssRes() "
					"logical error1" << std::endl);
		}
		status = sliding;
	}
	if (status == sliding) {
		// std::cerr << "v*current_velocity: " << v*current_velocity << std::endl;
		if (v*current_velocity < 0.) {
			if (((transition_type != from_sticked_to_sliding) &&
				(transition_type != from_sticking_to_sliding)) &&
				((d2DNorm(v-current_velocity) < d2DNorm(previous_switch_v)) ||
					(first_switch == true))) {
				// std::cerr << "XXXX" << std::endl;
				first_switch = false;
				status = sticking;
				transition_type = from_sliding_to_sticking;
				previous_switch_v = vel_ratio*(v-current_velocity);
				saved_sliding_velocity = v;
				saved_sliding_friction = f;
			}
		}
 	}

	switch (status) {
		case sticking: {
			//switch to sticking: null velocity at the end of time step
			current_friction_force = f;
			// std::cerr << "sticking"  << std::endl;
			WorkVec.IncCoef(startdof+1, v.x[0]);
			WorkVec.IncCoef(startdof+2, v.x[1]);
			break;
		}
		case sliding: {
			// std::cerr << "sliding"  << std::endl;
			doublereal vm = d2DNorm(v);
			//still sliding
			switch (transition_type) {
				case from_sticked_to_sliding: {
					// std::cerr << "from_sticked_to_sliding"  << std::endl;
					current_friction_force = fss(vm)*d2DDirection(f)+sigma2*v;
					break;
				}
				case from_sticking_to_sliding: {
					// std::cerr << "from_sticking_to_sliding"  << std::endl;
					current_friction_force = fss(vm) * d2DDirection(saved_sliding_friction) + sigma2 * v;
					break;
				}
				default: {
					if (vm >= vel_tolerance) {
						if (Dot(v, current_velocity) > 0.) {
							// std::cerr << "xx1"  << std::endl;
							// std::cerr << "v: "  << v << std::endl;
							// std::cerr << "curr v: " << current_velocity << std::endl;
							current_friction_force = fss(vm)*d2DDirection(v)+sigma2*v;
							use_sliding_v = true;
						} else {
							// std::cerr << "xx2"  << std::endl;
							// std::cerr << "v: "  << v << std::endl;
							// std::cerr << "curr v: " << current_velocity << std::endl;
							current_friction_force = fss(vm)*d2DDirection(f)+sigma2*v;
						}
					} else {
						// std::cerr << "xx3"  << std::endl;
						// std::cerr << "v: "  << v << std::endl;
						// std::cerr << "curr v: " << current_velocity << std::endl;
						//limit the force value while taking the sticking force d2DDirection
						current_friction_force = fss(vm)*d2DDirection(f)+sigma2*v;
					}
					if (vm < d2DNorm(current_velocity) && !first_iter) {
						// std::cerr << "xx4"  << std::endl;
						// std::cerr << "v: "  << v << std::endl;
						// std::cerr << "curr v: " << current_velocity << std::endl;
						current_velocity = v;
					}
					break;
				}
			}
			//save friction force value in the (algebric) state
			// std::cerr << "current_friction_force " << current_friction_force << std::endl;
			// std::cerr << "f " << f << std::endl;
			// std::cerr << "v " << v << std::endl;
			WorkVec.IncCoef(startdof+1, f.x[0] - current_friction_force.x[0]);
			WorkVec.IncCoef(startdof+2, f.x[1] - current_friction_force.x[1]);
			break;
		}
		case sticked: {
			// std::cerr << "sticked"  << std::endl;
			current_friction_force = f;
			WorkVec.IncCoef(startdof+1, v.x[0]);
			WorkVec.IncCoef(startdof+2, v.x[1]);
			break;
		}
		default: {
			silent_cerr("DiscreteCoulombFriction2D::AssRes() "
				"logical error" << std::endl);
		}
	}
	//update status
	first_iter = false;
	if (transition_type != null) {
		throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
	}
//	current_velocity = v;
};

void DiscreteCoulombFriction2D::AssJac(
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
	const ExpandableMatrix& dv) const {
	doublereal vm = d2DNorm(v);
	switch (status) {
		case sticking:
		case sticked: {
			//null velocity at the end of time step
			dv.SubFrom(WorkMat,startdof+1);
			dfc.ReDim(2, 1);
			dfc.SetBlockDim(1, 2);
			dfc.SetBlockIdx(1, startdof+1);
			dfc.Set(1., 1, 1, 1);
			dfc.Set(1., 2, 1, 2);
			break;
		}
		case sliding: {
			//still sliding
			//save friction force value in the (algebric) state
					// if (vm > 0.) {
					// 	if (Dot(v, current_velocity) > 0.) {
					// 		current_friction_force = fss(vm)*d2DDirection(v)+sigma2*v;
					// 	} else {
					// 		current_friction_force = fss(vm)*d2DDirection(f)+sigma2*v;
					// 	}
					// } else {
					// 	//limit the force value while taking the sticking force d2DDirection
					// 	current_friction_force = fss(vm)*d2DDirection(f)+sigma2*v;
					// }
					// if (vm < d2DNorm(current_velocity) && !first_iter) {
					// 	current_velocity = v;
					// }

			WorkMat.IncCoef(startdof+1,startdof+1,-1);
			WorkMat.IncCoef(startdof+2,startdof+2,-1);
			doublereal fssd = fss.ComputeDiff(vm);
			d2D dir = d2DDirection(current_friction_force-sigma2*v);
			d2D vm_d; d2DNorm_d(v, vm_d);

// 			if (use_sliding_v) {
// 				dfc.ReDim(2, 2);
// 				dfc.SetBlockDim(2, 2);
// 
// 				d2DDirection_d(v, Direction_d);
// 				Direction_d.Link(1, &dv);
// 				
// 				doublereal fs = fss(vm);
// 				dfc.Set(fs, 1, 2, 1);
// 				dfc.Set(fs, 2, 2, 2);
// 				dfc.Link(2, &Direction_d);
// 			} else {
// 				dfc.ReDim(2, 1);
// 			}
// 			dfc.SetBlockDim(1, 2);
// 			dfc.Set(fssd * dir.x[0] * vm_d.x[0] + sigma2, 1, 1, 1);
// 			dfc.Set(fssd * dir.x[0] * vm_d.x[1]         , 1, 1, 2);
// 			dfc.Set(fssd * dir.x[1] * vm_d.x[0]         , 2, 1, 1);
// 			dfc.Set(fssd * dir.x[1] * vm_d.x[1] + sigma2, 2, 1, 2);
// 			dfc.Link(1, &dv);
// 			dfc.Add(WorkMat, startdof+1, 1.);

			dfc.ReDim(2, 1);
			dfc.SetBlockDim(1, 2);
			dfc.Set(fssd * dir.x[0] * vm_d.x[0] + sigma2, 1, 1, 1);
			dfc.Set(fssd * dir.x[0] * vm_d.x[1]         , 1, 1, 2);
			dfc.Set(fssd * dir.x[1] * vm_d.x[0]         , 2, 1, 1);
			dfc.Set(fssd * dir.x[1] * vm_d.x[1] + sigma2, 2, 1, 2);
			dfc.Link(1, &dv);
			dfc.AddTo(WorkMat, startdof+1, 1.);
	
			// d2D diff = fss.ComputeDiff(vm)*d2DDirection(current_friction_force-sigma2*v)+sigma2*d2D({1., 1.});
			// dv.Add(WorkMat,startdof+1, diff.x[0]);
			// dv.Add(WorkMat,startdof+2, diff.x[1]);
			// dfc.ReDim(2, 1);
			// dfc.SetBlockDim(1, 2);
			// d2D diff2 = fss.ComputeDiff(vm)*d2DDirection(current_friction_force-sigma2*v)+sigma2*d2D({1., 1.});
			// dfc.Set(diff2.x[0], 1, 1, 1);
			// dfc.Set(diff2.x[1], 2, 1, 2);
			// dfc.Link(1, &dv);
			break;
		}
		default: {
			silent_cerr("DiscreteCoulombFriction2D::AssJac() "
				"logical error" << std::endl);
		}
	}
};

const OutputHandler::Dimensions
DiscreteCoulombFriction2D::GetEquationDimension(integer index) const {
	// DOF == 1

	OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
	case 1:
		dimension = OutputHandler::Dimensions::UnknownDimension;
		break;
	}

	return dimension;
}



//------------------------
d2D SimpleShapeCoefficient2D::Sh_c(void) const {
	return shc;
}
	
d2D SimpleShapeCoefficient2D::Sh_c(
	const d2D f,
	const doublereal F,
	const d2D v) {
	shc = f;
	return shc;
};

void SimpleShapeCoefficient2D::dSh_c(
	ExpandableMatrix& dShc,
	const d2D f,
	const doublereal F,
	const d2D v,
	const ExpandableMatrix& dfc,
	const ExpandableRowVector& dF,
	const ExpandableMatrix& dv) const {
		dShc.ReDim(2, 1);
		dShc.SetBlockDim(1, 2);
		dShc.Set(1., 1, 1, 1);
		dShc.Set(0., 1, 1, 2);
		dShc.Set(0., 2, 1, 1);
		dShc.Set(1., 2, 1, 2);
		dShc.Link(1, &dfc);
};



//---------------------------------------

BasicFriction2D *const ParseFriction2D(MBDynParser& HP,
	DataManager * pDM) 
{
   const char* sKeyWords[] = { 
      "modlugre" "2D",
      "discrete" "coulomb" "2D",
      NULL
   };
	enum KeyWords { 
	     MODLUGRE2D = 0,
	     DISCRETECOULOMB,
	     LASTKEYWORD
	};
	/* token corrente */
	KeyWords FuncType;
	
	KeyTable K(HP, sKeyWords);
	
	FuncType = KeyWords(HP.IsKeyWord());
	switch (FuncType) {
	case MODLUGRE2D: {
		doublereal sigma0 = HP.GetReal();
		doublereal sigma1 = HP.GetReal();
		doublereal sigma2 = HP.GetReal();
		doublereal kappa = HP.GetReal();
		const BasicScalarFunction*const sf =
			ParseScalarFunction(HP, pDM);
		return new ModLugreFriction2D(sigma0, sigma1, sigma2, kappa, sf);
		break;
	}
	case DISCRETECOULOMB: {
		const BasicScalarFunction*const sf =
			ParseScalarFunction(HP, pDM);
		doublereal sigma2 = 0.;
		doublereal vel_ratio = 0.8;
		doublereal vel_tolerance = 1.E-6;
		if (HP.IsKeyWord("sigma2")) {
			sigma2 = HP.GetReal();
		}
		if (HP.IsKeyWord("velocity" "ratio")) {
			vel_ratio = HP.GetReal();
		}
		if (HP.IsKeyWord("velocity" "tolerance")) {
			vel_tolerance = HP.GetReal();
		}
		return new DiscreteCoulombFriction2D(sf,sigma2, vel_ratio, vel_tolerance);
		break;
	}
	default: {
		silent_cerr("ParseFriction2D(): unrecognized friction type "
				"at line " << HP.GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		break;
	}
	}
	return 0;
};

BasicShapeCoefficient2D *const ParseShapeCoefficient2D(MBDynParser& HP) {
	const char* sKeyWords[] = {
		"simple" "2D",
 		NULL
	};
	enum KeyWords { 
		SIMPLE2D = 0,
		LASTKEYWORD
	};
	/* token corrente */
	KeyWords FuncType;
	
	KeyTable K(HP, sKeyWords);
	
	FuncType = KeyWords(HP.IsKeyWord());
	switch (FuncType) {
	case SIMPLE2D: {
		return new SimpleShapeCoefficient2D();
		break;
	}
	default: {
		silent_cerr("ParseShapeCoefficient2D(): "
			"unrecognized shape coefficient type "
			"at line " << HP.GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		break;
	}
	}
	return 0;
};
