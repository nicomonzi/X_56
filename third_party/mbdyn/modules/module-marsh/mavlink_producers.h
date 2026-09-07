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

#ifndef MODULE_MARSH_PRODUCERS_H
#define MODULE_MARSH_PRODUCERS_H

#include <cstring>
#include <string>
#include <vector>

/* pulls in dataman.h, which instruments.h and scalarvalue.h rely upon */
#include "mavlink_pubsub.h"

#include "instruments.h"
#include "scalarvalue.h"

/*
 * Where an outbound message takes the aircraft state from.
 *
 * Everything is returned in MAVLink conventions: body axes with x forward,
 * y right and z down, SI units, and NED for the "world" components.
 * Quantities a given source cannot provide are reported as unavailable,
 * and are left at zero unless the user overrides the corresponding field
 * (see MavlinkFieldDesc below).
 */
class MarshMotionSource {
public:
	virtual ~MarshMotionSource(void);

	/* node the state comes from, if any, so that the owning element can
	 * report it via GetConnectedNodes() */
	virtual const StructNode *pGetNode(void) const = 0;

	/* orientation of the aircraft frame */
	virtual Mat3x3 GetR(void) const = 0;
	/* the same orientation, as Euler parameters */
	virtual void GetQuaternion(doublereal& e0, Vec3& e) const;
	/* angular velocity, body axes */
	virtual Vec3 GetBodyRate(void) const = 0;
	/* specific force, namely the acceleration an accelerometer would
	 * measure (gravity included), body axes */
	virtual Vec3 GetSpecificForce(const Vec3& gravity) const;
	/* ground velocity, NED components */
	virtual Vec3 GetVelNED(void) const;

	/* Euler angles, only for sources that define them as flight
	 * mechanics does; the quaternion from GetR() is always available */
	virtual bool GetEulerAngles(Vec3& rpy) const {
		return false;
	};
	/* latitude and longitude in radians, altitude in m */
	virtual bool GetGeodetic(doublereal& lat, doublereal& lon,
		doublereal& alt) const
	{
		return false;
	};
};

/*
 * State of a structural node, with an optional offset and rotation; the
 * rotation doubles as the conversion from the axes of the node to MAVLink
 * body axes. Same pattern as modules/module-imu.
 */
class NodeMotionSource : public MarshMotionSource {
private:
	const StructNode *m_pNode;
	Vec3 m_tilde_f;
	Mat3x3 m_tilde_Rh;

public:
	NodeMotionSource(DataManager *pDM, MBDynParser& HP);

	const StructNode *pGetNode(void) const override;
	Mat3x3 GetR(void) const override;
	Vec3 GetBodyRate(void) const override;
	Vec3 GetSpecificForce(const Vec3& gravity) const override;
	Vec3 GetVelNED(void) const override;
};

/*
 * State of an "aircraft instruments" element: it already resolves the
 * aircraft frame, the flight mechanics angles and the geodetic position,
 * so it fills in the whole message.
 */
class InstrumentsMotionSource : public MarshMotionSource {
private:
	const AircraftInstruments *m_pInstruments;
	const StructNode *m_pNode;

	/* private data indices, resolved once while parsing */
	unsigned m_idxE0, m_idxE1, m_idxE2, m_idxE3;
	unsigned m_idxBank, m_idxAttitude, m_idxHeading;
	unsigned m_idxRollRate, m_idxPitchRate, m_idxYawRate;
	unsigned m_idxLatitude, m_idxLongitude, m_idxAltitude;
	unsigned m_idxVNorth, m_idxVEast, m_idxVDown;

	unsigned iGetIdx(const char *s) const;

public:
	InstrumentsMotionSource(DataManager *pDM, MBDynParser& HP);

	const StructNode *pGetNode(void) const override;
	Mat3x3 GetR(void) const override;
	void GetQuaternion(doublereal& e0, Vec3& e) const override;
	Vec3 GetBodyRate(void) const override;
	Vec3 GetVelNED(void) const override;
	bool GetEulerAngles(Vec3& rpy) const override;
	bool GetGeodetic(doublereal& lat, doublereal& lon,
		doublereal& alt) const override;
};

/*
 * Reads the source of an outbound message, either
 *     aircraft instruments , <label>
 * or
 *     node , <label> [ , position , ... ] [ , orientation , ... ]
 * Returns 0, without consuming anything, if neither keyword is there: a
 * message may well be assembled out of drive callers alone.
 */
extern MarshMotionSource *ReadMarshMotionSource(DataManager *pDM,
	MBDynParser& HP);

/*
 * Per-field overrides: any (float) field of a message can be driven by an
 * arbitrary drive caller or node dof, either to complete what the source
 * leaves out, or to inject something the aircraft dynamics does not
 * describe (special effects, faults, ...). An override always wins over
 * the value computed from the source.
 */
template <class M>
struct MavlinkFieldDesc {
	const char *name;
	float M::*pMember;
};

template <class M>
struct MavlinkFieldValue {
	float M::*pMember;
	ScalarValue *pValue;
};

/* reads zero or more "field , <name> , <value>" clauses */
template <class M>
void
ReadMavlinkFieldOverrides(DataManager *pDM, MBDynParser& HP,
	const MavlinkFieldDesc<M> *desc,
	std::vector<MavlinkFieldValue<M> >& fields)
{
	while (HP.IsKeyWord("field")) {
		const std::string name = HP.GetStringWithDelims();

		float M::*pMember = 0;
		for (unsigned i = 0; desc[i].name != 0; i++) {
			if (name == desc[i].name) {
				pMember = desc[i].pMember;
				break;
			}
		}

		if (pMember == 0) {
			silent_cerr("module-marsh: unknown message field \""
				<< name << "\" at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		MavlinkFieldValue<M> f;
		f.pMember = pMember;
		f.pValue = ReadScalarValue(pDM, HP);
		fields.push_back(f);
	}
}

template <class M>
void
ApplyMavlinkFieldOverrides(M& msg,
	const std::vector<MavlinkFieldValue<M> >& fields)
{
	for (typename std::vector<MavlinkFieldValue<M> >::const_iterator
		i = fields.begin(); i != fields.end(); ++i)
	{
		msg.*(i->pMember) = i->pValue->dGetValue();
	}
}

class SimStateProducer : public MavlinkProducer {
private:
	MarshMotionSource *m_pSource;
	std::vector<MavlinkFieldValue<mavlink_sim_state_t> > m_fields;

public:
	SimStateProducer(DataManager *pDM, MBDynParser& HP);
	virtual ~SimStateProducer(void);

	const StructNode *pGetNode(void) const override;
	void Pack(mavlink_message_t& msg, uint32_t time_boot_ms,
		uint8_t sysid, uint8_t compid, const Vec3& gravity) const override;
};

class MotionCueExtraProducer : public MavlinkProducer {
private:
	MarshMotionSource *m_pSource;
	std::vector<MavlinkFieldValue<mavlink_motion_cue_extra_t> > m_fields;

public:
	MotionCueExtraProducer(DataManager *pDM, MBDynParser& HP);
	virtual ~MotionCueExtraProducer(void);

	const StructNode *pGetNode(void) const override;
	void Pack(mavlink_message_t& msg, uint32_t time_boot_ms,
		uint8_t sysid, uint8_t compid, const Vec3& gravity) const override;
};

#endif /* MODULE_MARSH_PRODUCERS_H */
