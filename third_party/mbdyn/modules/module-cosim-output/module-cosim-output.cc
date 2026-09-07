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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

//#include <cstdint>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <fstream>

#include "userelem.h"
#include "socketstream_out_elem.h"
#include "socketstreamdrive.h"
#include "bufmod.h"
//#include "drive_.h"
#include "stepsol.h"


// StreamContentCosim - begin

class StreamContentCosim : public StreamContent {
protected:
	const DataManager *m_pDM;
	const StructNode *m_pNode;
	doublereal m_h;
	doublereal m_dRho;

public:
	StreamContentCosim(const DataManager *pDM, const StructNode *pNode,
		doublereal h, doublereal rho,
		StreamContent::Modifier *pMod /* not needed */ );
	virtual ~StreamContentCosim(void);

	void Prepare(void);
	unsigned GetNumChannels(void) const;
};

// #include "geomdata.h"

StreamContentCosim::StreamContentCosim(const DataManager *pDM,
	const StructNode * pNode,
	doublereal h, doublereal rho,
	StreamContent::Modifier *pMod)
: StreamContent(0, pMod), m_pDM(pDM), m_pNode(pNode), m_h(h), m_dRho(rho)
{
	unsigned int size = sizeof(doublereal)*6;

	ASSERT(h > 0.);
	ASSERT((rho >= 0.) and (rho <= 1.));

	buf.resize(size);
	memset(&buf[0], 0, size);
	m_pMod->Set(size, &buf[0]);
}

StreamContentCosim::~StreamContentCosim(void)
{
	NO_OP;
}

void
StreamContentCosim::Prepare(void)
{
	doublereal h = m_h;
	if (h == -1.) {
		// NOTE: this is the current time step, not the next one;
		// we are assuming the problem is integrated with constant time step
		h = m_pDM->pGetDrvHdl()->dGetTimeStep();
	}

	integer iStep = m_pDM->pGetDrvHdl()->iGetStep();
	if (iStep == 0) {
		// first step: nothing is known about the previous step
		// use Crank Nicolson, extrapolating constant velocity (explicit Euler)

		doublereal *dbuf = (doublereal *)&buf[0];

		// predicted position
		const Vec3& X(m_pNode->GetXCurr());
		const Vec3& V(m_pNode->GetVCurr());

		dbuf[0] = X(1) + V(1)*h;
		dbuf[1] = X(2) + V(2)*h;
		dbuf[2] = X(3) + V(3)*h;

		// predicted orientation
		const Mat3x3& R = m_pNode->GetRCurr();
		const Vec3& Omega = m_pNode->GetWCurr();
		Vec3 g_p1(Omega*h);
		Mat3x3 RDelta = Mat3x3(CGR_Rot::MatR, g_p1);
		Vec3 Theta(RotManip::VecRot(RDelta*R));

		dbuf[3 + 0] = Theta(1);
		dbuf[3 + 1] = Theta(2);
		dbuf[3 + 2] = Theta(3);

	} else if (iStep > 0) {
		// for subsequent steps, the state at the previous steps is available

		// Prediction (using fixed step h and linear two-step method with tunable algorithmic dissipation
		doublereal dAlpha = 1.;
		doublereal dDen = 2.*(1. + dAlpha) - (1. - m_dRho)*(1. - m_dRho);
		doublereal dBeta = dAlpha*((1. - m_dRho)*(1. - m_dRho)*(2. + dAlpha)
			+2.*(2.*m_dRho - 1.)*(1. + dAlpha))/dDen;
		doublereal dDelta = .5*dAlpha*dAlpha*(1. - m_dRho)*(1. - m_dRho)/dDen;

		doublereal m_mp[2], m_np[2], m_a[StepNIntegrator::IDX_A2 + 1], m_b[StepNIntegrator::IDX_B2 + 1];

		m_mp[0] = -6.*dAlpha*dAlpha*(1. + dAlpha);
		m_mp[1] = -m_mp[0];
		m_np[0] = 1. + 4.*dAlpha + 3.*dAlpha*dAlpha;
		m_np[1] = dAlpha*(2. + 3.*dAlpha);

		m_a[StepNIntegrator::IDX_A1] = 1. - dBeta;
		m_a[StepNIntegrator::IDX_A2] = dBeta;
		m_b[StepNIntegrator::IDX_B0] = h*(dDelta/dAlpha + dAlpha/2.);
		m_b[StepNIntegrator::IDX_B1] = h*(dBeta/2. + dAlpha/2. - dDelta/dAlpha*(1. + dAlpha));
		m_b[StepNIntegrator::IDX_B2] = h*(dBeta/2. + dDelta);

		m_mp[0] /= h;
		m_mp[1] /= h;

		doublereal *dbuf = (doublereal *)&buf[0];

		// predicted position
		const Vec3& X = m_pNode->GetXCurr();
		const Vec3& X_m1 = m_pNode->GetXPrev();
		const Vec3& V = m_pNode->GetVCurr();
		const Vec3& V_m1 = m_pNode->GetVPrev();

		for (integer idx = 0; idx < 3; ++idx) {
			// predicted velocity component (derivative of cubic extrapolation on position)
			doublereal dx_p1 = m_mp[0]*X(idx + 1) + m_mp[1]*X_m1(idx + 1) + m_np[0]*V(idx + 1) + m_np[1]*V_m1(idx + 1);

			// predicted position (two-step linear multistep method with algorithmic dissipation)
			dbuf[idx] = m_a[StepNIntegrator::IDX_A1]*X(idx + 1) + m_a[StepNIntegrator::IDX_A2]*X_m1(idx + 1)
				+ m_b[StepNIntegrator::IDX_B0]*dx_p1
				+ m_b[StepNIntegrator::IDX_B1]*V(idx + 1) + m_b[StepNIntegrator::IDX_B2]*V_m1(idx + 1);
		}

		// predicted orientation
		const Mat3x3& R = m_pNode->GetRCurr();
		const Mat3x3& R_m1 = m_pNode->GetRPrev();
		const Vec3& Omega = m_pNode->GetWCurr();
		const Vec3& Omega_m1 = m_pNode->GetWPrev();

		// backward relative rotation between current and previous step
		Mat3x3 RDelta(R_m1.MulMT(R));
		// assuming g == 0
		Vec3 g_m1(CGR_Rot::Param, RDelta);
		Vec3 gP_m1(Mat3x3(CGR_Rot::MatGm1, g_m1)*Omega_m1);
		const Vec3& gP(Omega);

		Vec3 g_p1;
		for (integer idx = 0; idx < 3; ++idx) {
			// predicted orientation parameter derivative component (derivative of cubic extrapolation on position)
			doublereal dg_p1 = /* m_mp[0]*g(idx + 1) + */ m_mp[1]*g_m1(idx + 1) + m_np[0]*gP(idx + 1) + m_np[1]*gP_m1(idx + 1);

			// predicted orientation parameter (two-step linear multistep method with algorithmic dissipation)
			g_p1(idx + 1) = /* m_a[StepNIntegrator::IDX_A1]*g(idx + 1) + */ m_a[StepNIntegrator::IDX_A2]*g_m1(idx + 1)
				+ m_b[StepNIntegrator::IDX_B0]*dg_p1
				+ m_b[StepNIntegrator::IDX_B1]*gP(idx + 1) + m_b[StepNIntegrator::IDX_B2]*gP_m1(idx + 1);
		}

		// forward relative rotation between current and next step
		RDelta = Mat3x3(CGR_Rot::MatR, g_p1);
		// absolute Euler vector of predicted orientation
		Vec3 Theta(RotManip::VecRot(RDelta*R));

		dbuf[3 + 0] = Theta(1);
		dbuf[3 + 1] = Theta(2);
		dbuf[3 + 2] = Theta(3);
	}

	m_pMod->Modify();
}

unsigned
StreamContentCosim::GetNumChannels(void) const
{
	return buf.size()/sizeof(doublereal);
}

// StreamContentCosim - end


// Read stream output element content type for cosimulation
struct CosimStreamOutputReader : public StreamOutputContentTypeReader {
	virtual StreamContent* Read(DataManager* pDM, MBDynParser& HP);
};

StreamContent*
CosimStreamOutputReader::Read(DataManager* pDM, MBDynParser& HP)
{
	StreamContent* pSC(0);

	/*
	   Syntax:

	   cosimulation,
	       <node_label>,
	       time step, { <time_step> | from solver }
	       [ , rho, <rho> ]

	   when "from solver" is used, the time step is obtained from the DataManager;
	   NOTE: this corresponds to the *current* time step, not the one
           that will be used by the solver for the actual prediction.
	   Otherwise, the time step must be greater than zero;
	   it is the user's responsibility to make sure it matches the time step used
	   by the solver

	   If rho is not given, it defaults to zero (second-order BDF); otherwise, it
	   must be between 0 and 1

	 */

	const StructNode* pNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	// -1 means that it will be taken from the solver!
	doublereal h(-1.);
	if (HP.IsKeyWord("time" "step")) {
		if (HP.IsKeyWord("from" "solver")) {
			silent_cout("CosimStreamOutput: the time step will be taken from the simulation at line " << HP.GetLineData() << "; "
				"note that the time step of the current step will be used, not that of the subsequent step; "
				"thus, this only works when the time step is constant." << std::endl);

		} else {
			try {
				h = HP.GetReal(0., HighParser::range_gt<doublereal>(0.));

			} catch (HighParser::ErrValueOutOfRange<doublereal>& e) {
				silent_cerr("error: invalid reference time step" << e.Get() << " (must be positive) [" << e.what() << "] for CosimStreamOutput at line " << HP.GetLineData() << std::endl);
				throw e;
			}

			silent_cout("CosimStreamOutput: a constant time step h=" << h << " will be used at line " << HP.GetLineData() << "; "
				"make sure it matches the time step used in the simulation." << std::endl);
		}

	} else {
		silent_cerr("time step required for CosimStreamOutput at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal rho(0.);
	if (HP.IsKeyWord("bdf") || HP.IsKeyWord("bdf2")) {
		// BDF2, the default
		rho = 0.;

	} else if (HP.IsKeyWord("rho")) {
		try {
			rho = HP.GetReal(0., HighParser::range_ge_le<doublereal>(0., 1.));

		} catch (HighParser::ErrValueOutOfRange<doublereal>& e) {
			silent_cerr("error: invalid reference asymptotic spectral radius rho=" << e.Get() << " (must be greater than or equal to zero and lower than or equal to one) [" << e.what() << "] for CosimStreamOutput at line " << HP.GetLineData() << std::endl);
			throw e;
		}

	} else {
		silent_cout("Warning, no spectral radius provided; using rho=0 (BDF2) at line " << HP.GetLineData() << std::endl);
	}

	StreamContent::Modifier *pMod(0);

	SAFENEWWITHCONSTRUCTOR(pSC, StreamContentCosim, StreamContentCosim(pDM, pNode, h, rho, pMod));

	return pSC;
}

extern "C"
int module_init(const char *module_name, void *pdm, void *php){

#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif
	// Stream output element content type for cosimulation
	StreamOutputContentTypeReader *rf = new CosimStreamOutputReader;
	if (!SetStreamOutputContentType("cosimulation", rf)) {
		delete rf;
		return -1;
	}

	return 0;
}

