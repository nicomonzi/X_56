/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2026
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 */

#include "mbconfig.h"

#include <cmath>
#include <cstring>
#include <map>

#include "mavlink_producers.h"

/* wraps an angle to (-pi, pi], as MAVLink expects Euler angles */
static doublereal
WrapAngle(doublereal a)
{
	while (a > M_PI) {
		a -= 2*M_PI;
	}
	while (a <= -M_PI) {
		a += 2*M_PI;
	}

	return a;
}

/*
 * Rotation matrix out of Euler parameters; the inverse of
 * MatR2EulerParams(), which MBDyn does not provide.
 */
static Mat3x3
EulerParams2MatR(const doublereal& e0, const Vec3& e)
{
	return ::Eye3*(e0*e0 - e.Dot()) + e.Tens()*2. + Mat3x3(MatCross, e)*(2.*e0);
}

/* MarshMotionSource - begin */

MarshMotionSource::~MarshMotionSource(void)
{
	NO_OP;
}

void
MarshMotionSource::GetQuaternion(doublereal& e0, Vec3& e) const
{
	MatR2EulerParams(GetR(), e0, e);
}

Vec3
MarshMotionSource::GetSpecificForce(const Vec3& gravity) const
{
	/* MAVLink reports the acceleration as specific force, namely what an
	 * accelerometer measures, which includes gravity: a body at rest
	 * reports a vector whose magnitude is that of gravity, not a null
	 * one. This is the mirror image of the subtraction the MARSH
	 * reference node (mps-adapter) does when consuming these messages. */
	const StructNode *pNode = pGetNode();
	if (pNode == 0) {
		return ::Zero3;
	}

	return GetR().MulTV(pNode->GetXPPCurr() - gravity);
}

Vec3
MarshMotionSource::GetVelNED(void) const
{
	const StructNode *pNode = pGetNode();
	if (pNode == 0) {
		return ::Zero3;
	}

	/* MBDyn's flight mechanics "world" is North-West-Up, MAVLink's is
	 * North-East-Down */
	const Vec3& V(pNode->GetVCurr());

	return Vec3(V(1), -V(2), -V(3));
}

/* MarshMotionSource - end */

/* NodeMotionSource - begin */

NodeMotionSource::NodeMotionSource(DataManager *pDM, MBDynParser& HP)
: m_pNode(0),
m_tilde_f(::Zero3),
m_tilde_Rh(::Eye3)
{
	m_pNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
	if (!m_pNode->bComputeAccelerations()) {
		const_cast<StructNode *>(m_pNode)->ComputeAccelerations(true);
	}

	ReferenceFrame RF(m_pNode);
	if (HP.IsKeyWord("position")) {
		m_tilde_f = HP.GetPosRel(RF);
	}

	if (HP.IsKeyWord("orientation")) {
		m_tilde_Rh = HP.GetRotRel(RF);
	}
}

const StructNode *
NodeMotionSource::pGetNode(void) const
{
	return m_pNode;
}

Mat3x3
NodeMotionSource::GetR(void) const
{
	return m_pNode->GetRCurr()*m_tilde_Rh;
}

Vec3
NodeMotionSource::GetBodyRate(void) const
{
	return GetR().MulTV(m_pNode->GetWCurr());
}

Vec3
NodeMotionSource::GetSpecificForce(const Vec3& gravity) const
{
	const Mat3x3 R(GetR());
	Vec3 a(m_pNode->GetXPPCurr() - gravity);

	/* contribution of the offset, as in module-imu */
	const Vec3 overline_f(m_tilde_Rh.MulTV(m_tilde_f));
	if (!overline_f.IsNull()) {
		const Vec3 Omega(R.MulTV(m_pNode->GetWCurr()));
		const Vec3 OmegaP(R.MulTV(m_pNode->GetWPCurr()));

		return R.MulTV(a) + OmegaP.Cross(overline_f)
			+ Omega.Cross(Omega.Cross(overline_f));
	}

	return R.MulTV(a);
}

Vec3
NodeMotionSource::GetVelNED(void) const
{
	Vec3 V(m_pNode->GetVCurr());

	const Vec3 f(m_pNode->GetRCurr()*m_tilde_f);
	if (!f.IsNull()) {
		V += m_pNode->GetWCurr().Cross(f);
	}

	return Vec3(V(1), -V(2), -V(3));
}

/* NodeMotionSource - end */

/* InstrumentsMotionSource - begin */

InstrumentsMotionSource::InstrumentsMotionSource(DataManager *pDM, MBDynParser& HP)
: m_pInstruments(0),
m_pNode(0)
{
	const unsigned uLabel = HP.GetInt();

	const Elem *pElem = pDM->pFindElem(Elem::AERODYNAMIC, uLabel);
	if (pElem == 0) {
		silent_cerr("module-marsh: unable to find AerodynamicElement("
			<< uLabel << ") at line " << HP.GetLineData()
			<< "; note that it must be declared before this element"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_pInstruments = dynamic_cast<const AircraftInstruments *>(pElem);
	if (m_pInstruments == 0) {
		silent_cerr("module-marsh: AerodynamicElement(" << uLabel
			<< ") is not an \"aircraft instruments\" element at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* the element only connects the node representing the aircraft */
	std::vector<const Node *> nodes;
	m_pInstruments->GetConnectedNodes(nodes);
	if (nodes.size() == 1) {
		m_pNode = dynamic_cast<const StructNode *>(nodes[0]);
	}

	if (m_pNode == 0) {
		silent_cerr("module-marsh: unable to get the structural node of "
			"AircraftInstruments(" << uLabel << ") at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!m_pNode->bComputeAccelerations()) {
		const_cast<StructNode *>(m_pNode)->ComputeAccelerations(true);
	}

	m_idxE0 = iGetIdx("e0");
	m_idxE1 = iGetIdx("e1");
	m_idxE2 = iGetIdx("e2");
	m_idxE3 = iGetIdx("e3");

	m_idxBank = iGetIdx("bank");
	m_idxAttitude = iGetIdx("attitude");
	m_idxHeading = iGetIdx("heading");

	m_idxRollRate = iGetIdx("rollrate");
	m_idxPitchRate = iGetIdx("pitchrate");
	m_idxYawRate = iGetIdx("yawrate");

	m_idxLatitude = iGetIdx("latitude");
	m_idxLongitude = iGetIdx("longitude");
	m_idxAltitude = iGetIdx("altitude");

	m_idxVNorth = iGetIdx("vnorth");
	m_idxVEast = iGetIdx("veast");
	m_idxVDown = iGetIdx("vdown");
}

unsigned
InstrumentsMotionSource::iGetIdx(const char *s) const
{
	const unsigned idx = m_pInstruments->iGetPrivDataIdx(s);
	if (idx == 0) {
		silent_cerr("module-marsh: AircraftInstruments("
			<< m_pInstruments->GetLabel() << ") does not provide "
			"private data \"" << s << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return idx;
}

const StructNode *
InstrumentsMotionSource::pGetNode(void) const
{
	return m_pNode;
}

void
InstrumentsMotionSource::GetQuaternion(doublereal& e0, Vec3& e) const
{
	/* the element resolves the aircraft frame, including the orientation
	 * of the aircraft with respect to its node, and exposes it directly
	 * as Euler parameters */
	e0 = m_pInstruments->dGetPrivData(m_idxE0);
	e = Vec3(m_pInstruments->dGetPrivData(m_idxE1),
		m_pInstruments->dGetPrivData(m_idxE2),
		m_pInstruments->dGetPrivData(m_idxE3));
}

Mat3x3
InstrumentsMotionSource::GetR(void) const
{
	doublereal e0;
	Vec3 e;
	GetQuaternion(e0, e);

	return EulerParams2MatR(e0, e);
}

Vec3
InstrumentsMotionSource::GetBodyRate(void) const
{
	return Vec3(m_pInstruments->dGetPrivData(m_idxRollRate),
		m_pInstruments->dGetPrivData(m_idxPitchRate),
		m_pInstruments->dGetPrivData(m_idxYawRate));
}

Vec3
InstrumentsMotionSource::GetVelNED(void) const
{
	return Vec3(m_pInstruments->dGetPrivData(m_idxVNorth),
		m_pInstruments->dGetPrivData(m_idxVEast),
		m_pInstruments->dGetPrivData(m_idxVDown));
}

bool
InstrumentsMotionSource::GetEulerAngles(Vec3& rpy) const
{
	rpy = Vec3(m_pInstruments->dGetPrivData(m_idxBank),
		m_pInstruments->dGetPrivData(m_idxAttitude),
		/* the element reports the heading in [0, 2*pi) */
		WrapAngle(m_pInstruments->dGetPrivData(m_idxHeading)));

	return true;
}

bool
InstrumentsMotionSource::GetGeodetic(doublereal& lat, doublereal& lon,
	doublereal& alt) const
{
	lat = m_pInstruments->dGetPrivData(m_idxLatitude);
	lon = m_pInstruments->dGetPrivData(m_idxLongitude);
	alt = m_pInstruments->dGetPrivData(m_idxAltitude);

	return true;
}

/* InstrumentsMotionSource - end */

MarshMotionSource *
ReadMarshMotionSource(DataManager *pDM, MBDynParser& HP)
{
	if (HP.IsKeyWord("aircraft" "instruments")) {
		return new InstrumentsMotionSource(pDM, HP);
	}

	if (HP.IsKeyWord("node")) {
		return new NodeMotionSource(pDM, HP);
	}

	return 0;
}

/* SimStateProducer - begin */

static const MavlinkFieldDesc<mavlink_sim_state_t> SimStateFields[] = {
	{ "q1", &mavlink_sim_state_t::q1 },
	{ "q2", &mavlink_sim_state_t::q2 },
	{ "q3", &mavlink_sim_state_t::q3 },
	{ "q4", &mavlink_sim_state_t::q4 },
	{ "roll", &mavlink_sim_state_t::roll },
	{ "pitch", &mavlink_sim_state_t::pitch },
	{ "yaw", &mavlink_sim_state_t::yaw },
	{ "xacc", &mavlink_sim_state_t::xacc },
	{ "yacc", &mavlink_sim_state_t::yacc },
	{ "zacc", &mavlink_sim_state_t::zacc },
	{ "xgyro", &mavlink_sim_state_t::xgyro },
	{ "ygyro", &mavlink_sim_state_t::ygyro },
	{ "zgyro", &mavlink_sim_state_t::zgyro },
	{ "lat", &mavlink_sim_state_t::lat },
	{ "lon", &mavlink_sim_state_t::lon },
	{ "alt", &mavlink_sim_state_t::alt },
	{ "std_dev_horz", &mavlink_sim_state_t::std_dev_horz },
	{ "std_dev_vert", &mavlink_sim_state_t::std_dev_vert },
	{ "vn", &mavlink_sim_state_t::vn },
	{ "ve", &mavlink_sim_state_t::ve },
	{ "vd", &mavlink_sim_state_t::vd },

	{ 0, 0 }
};

SimStateProducer::SimStateProducer(DataManager *pDM, MBDynParser& HP)
: m_pSource(0)
{
	m_pSource = ReadMarshMotionSource(pDM, HP);
	ReadMavlinkFieldOverrides(pDM, HP, SimStateFields, m_fields);
}

SimStateProducer::~SimStateProducer(void)
{
	if (m_pSource != 0) {
		delete m_pSource;
	}
}

const StructNode *
SimStateProducer::pGetNode(void) const
{
	return m_pSource != 0 ? m_pSource->pGetNode() : 0;
}

void
SimStateProducer::Pack(mavlink_message_t& msg, uint32_t /* time_boot_ms */,
	uint8_t sysid, uint8_t compid, const Vec3& gravity) const
{
	mavlink_sim_state_t state;
	memset(&state, 0, sizeof(state));

	if (m_pSource != 0) {
		doublereal e0;
		Vec3 e;
		m_pSource->GetQuaternion(e0, e);
		state.q1 = e0;
		state.q2 = e(1);
		state.q3 = e(2);
		state.q4 = e(3);

		Vec3 rpy;
		if (m_pSource->GetEulerAngles(rpy)) {
			state.roll = rpy(1);
			state.pitch = rpy(2);
			state.yaw = rpy(3);
		}

		const Vec3 accel(m_pSource->GetSpecificForce(gravity));
		state.xacc = accel(1);
		state.yacc = accel(2);
		state.zacc = accel(3);

		const Vec3 rate(m_pSource->GetBodyRate());
		state.xgyro = rate(1);
		state.ygyro = rate(2);
		state.zgyro = rate(3);

		doublereal lat, lon, alt;
		if (m_pSource->GetGeodetic(lat, lon, alt)) {
			/* the element works in radians, MAVLink in degrees */
			lat *= 180./M_PI;
			lon *= 180./M_PI;

			state.lat = lat;
			state.lon = lon;
			state.alt = alt;

			state.lat_int = int32_t(std::round(lat*1e7));
			state.lon_int = int32_t(std::round(lon*1e7));
		}

		const Vec3 vned(m_pSource->GetVelNED());
		state.vn = vned(1);
		state.ve = vned(2);
		state.vd = vned(3);
	}

	ApplyMavlinkFieldOverrides(state, m_fields);

	/* keep the high precision geodetic fields consistent with any
	 * override of the low precision ones */
	if (!m_fields.empty()) {
		state.lat_int = int32_t(std::round(doublereal(state.lat)*1e7));
		state.lon_int = int32_t(std::round(doublereal(state.lon)*1e7));
	}

	mavlink_msg_sim_state_encode(sysid, compid, &msg, &state);
}

/* SimStateProducer - end */

/* MotionCueExtraProducer - begin */

static const MavlinkFieldDesc<mavlink_motion_cue_extra_t> MotionCueExtraFields[] = {
	{ "vel_roll", &mavlink_motion_cue_extra_t::vel_roll },
	{ "vel_pitch", &mavlink_motion_cue_extra_t::vel_pitch },
	{ "vel_yaw", &mavlink_motion_cue_extra_t::vel_yaw },
	{ "acc_x", &mavlink_motion_cue_extra_t::acc_x },
	{ "acc_y", &mavlink_motion_cue_extra_t::acc_y },
	{ "acc_z", &mavlink_motion_cue_extra_t::acc_z },

	{ 0, 0 }
};

MotionCueExtraProducer::MotionCueExtraProducer(DataManager *pDM, MBDynParser& HP)
: m_pSource(0)
{
	m_pSource = ReadMarshMotionSource(pDM, HP);
	ReadMavlinkFieldOverrides(pDM, HP, MotionCueExtraFields, m_fields);
}

MotionCueExtraProducer::~MotionCueExtraProducer(void)
{
	if (m_pSource != 0) {
		delete m_pSource;
	}
}

const StructNode *
MotionCueExtraProducer::pGetNode(void) const
{
	return m_pSource != 0 ? m_pSource->pGetNode() : 0;
}

void
MotionCueExtraProducer::Pack(mavlink_message_t& msg, uint32_t time_boot_ms,
	uint8_t sysid, uint8_t compid, const Vec3& gravity) const
{
	mavlink_motion_cue_extra_t extra;
	memset(&extra, 0, sizeof(extra));

	extra.time_boot_ms = time_boot_ms;

	if (m_pSource != 0) {
		const Vec3 rate(m_pSource->GetBodyRate());
		extra.vel_roll = rate(1);
		extra.vel_pitch = rate(2);
		extra.vel_yaw = rate(3);

		const Vec3 accel(m_pSource->GetSpecificForce(gravity));
		extra.acc_x = accel(1);
		extra.acc_y = accel(2);
		extra.acc_z = accel(3);
	}

	ApplyMavlinkFieldOverrides(extra, m_fields);

	mavlink_msg_motion_cue_extra_encode(sysid, compid, &msg, &extra);
}

/* MotionCueExtraProducer - end */

/* registry - begin */

typedef MavlinkProducer *(*ProducerFactory)(DataManager *, MBDynParser&);

template <class T>
static MavlinkProducer *
ReadProducer(DataManager *pDM, MBDynParser& HP)
{
	return new T(pDM, HP);
}

static const std::map<std::string, ProducerFactory> ProducerRegistry = {
	{ "SIM_STATE", &ReadProducer<SimStateProducer> },
	{ "MOTION_CUE_EXTRA", &ReadProducer<MotionCueExtraProducer> },
};

MavlinkProducer *
ReadMarshProducer(const std::string& name, DataManager *pDM, MBDynParser& HP)
{
	std::map<std::string, ProducerFactory>::const_iterator i = ProducerRegistry.find(name);
	if (i == ProducerRegistry.end()) {
		silent_cerr("module-marsh: unknown outbound message \"" << name
			<< "\" in \"send\" clause at line " << HP.GetLineData()
			<< std::endl);
		return 0;
	}

	return (*i->second)(pDM, HP);
}

/* registry - end */
