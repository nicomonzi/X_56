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

#include <cstring>
#include <map>

#include "mavlink_consumers.h"

/* ManualControlConsumer - begin */

ManualControlConsumer::ManualControlConsumer(DataManager * /* pDM */, MBDynParser& HP)
: m_x(0), m_y(0), m_z(0), m_r(0), m_buttons(0)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"subscribe: \"MANUAL_CONTROL\";					\n"
"	exposes the decoded axes/buttons as private data:		\n"
"	x, y, z, r, buttons						\n"
			<< std::endl);
	}
}

uint32_t
ManualControlConsumer::GetMsgId(void) const
{
	return MAVLINK_MSG_ID_MANUAL_CONTROL;
}

void
ManualControlConsumer::Decode(const mavlink_message_t& msg)
{
	mavlink_manual_control_t mc;
	mavlink_msg_manual_control_decode(&msg, &mc);

	m_x = mc.x;
	m_y = mc.y;
	m_z = mc.z;
	m_r = mc.r;
	m_buttons = mc.buttons;
}

unsigned int
ManualControlConsumer::iGetNumPrivData(void) const
{
	return 5;
}

unsigned int
ManualControlConsumer::iGetPrivDataIdx(const char *s) const
{
	if (strcmp(s, "MANUAL_CONTROL.x") == 0) { return 1; }
	if (strcmp(s, "MANUAL_CONTROL.y") == 0) { return 2; }
	if (strcmp(s, "MANUAL_CONTROL.z") == 0) { return 3; }
	if (strcmp(s, "MANUAL_CONTROL.r") == 0) { return 4; }
	if (strcmp(s, "MANUAL_CONTROL.buttons") == 0) { return 5; }

	return 0;
}

doublereal
ManualControlConsumer::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1: return doublereal(m_x);
	case 2: return doublereal(m_y);
	case 3: return doublereal(m_z);
	case 4: return doublereal(m_r);
	case 5: return doublereal(m_buttons);
	}

	return 0.;
}

/* ManualControlConsumer - end */

/* registry - begin */

typedef MavlinkConsumer *(*ConsumerFactory)(DataManager *, MBDynParser&);

template <class T>
static MavlinkConsumer *
ReadConsumer(DataManager *pDM, MBDynParser& HP)
{
	return new T(pDM, HP);
}

static const std::map<std::string, ConsumerFactory> ConsumerRegistry = {
	{ "MANUAL_CONTROL", &ReadConsumer<ManualControlConsumer> },
};

MavlinkConsumer *
ReadMarshConsumer(const std::string& name, DataManager *pDM, MBDynParser& HP)
{
	std::map<std::string, ConsumerFactory>::const_iterator i = ConsumerRegistry.find(name);
	if (i == ConsumerRegistry.end()) {
		silent_cerr("module-marsh: unknown inbound message \"" << name
			<< "\" in \"subscribe\" clause at line " << HP.GetLineData()
			<< std::endl);
		return 0;
	}

	return (*i->second)(pDM, HP);
}

/* registry - end */
