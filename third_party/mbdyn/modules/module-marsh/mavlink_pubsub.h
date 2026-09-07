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

#ifndef MODULE_MARSH_PUBSUB_H
#define MODULE_MARSH_PUBSUB_H

#include <cstdint>
#include <string>

#include "dataman.h"
#include "userelem.h"

#include "mavlink/marsh/mavlink.h"

/*
 * Something that knows how to fill in one outbound MAVLink message from
 * live MBDyn state. Registered by name (e.g. "SIM_STATE") so a `send`
 * clause in the .mbd deck can select it.
 */
class MavlinkProducer {
public:
	virtual ~MavlinkProducer(void) { };

	/* fills in msg; does not send it */
	virtual void Pack(mavlink_message_t& msg, uint32_t time_boot_ms,
		uint8_t sysid, uint8_t compid, const Vec3& gravity) const = 0;

	/* node this producer reads from, if any, so the owning element can
	 * report it via GetConnectedNodes() */
	virtual const StructNode *pGetNode(void) const {
		return 0;
	};
};

/*
 * Something that knows how to decode one inbound MAVLink message and,
 * optionally, expose the decoded values as named private data. Registered
 * by name (e.g. "MANUAL_CONTROL") so a `subscribe` clause in the .mbd deck
 * can select it.
 */
class MavlinkConsumer {
public:
	virtual ~MavlinkConsumer(void) { };

	virtual uint32_t GetMsgId(void) const = 0;
	virtual void Decode(const mavlink_message_t& msg) = 0;

	virtual unsigned int iGetNumPrivData(void) const {
		return 0;
	};
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return 0;
	};
	virtual doublereal dGetPrivData(unsigned int i) const {
		return 0.;
	};
};

/*
 * Registry lookup: build a Producer/Consumer instance by MAVLink message
 * name, consuming any message-specific arguments from HP. Returns 0 (after
 * printing an error) if name is not a known message.
 *
 * Adding a new message type is: write the subclass, add one line to the
 * table in the corresponding .cc file.
 */
extern MavlinkProducer *ReadMarshProducer(const std::string& name,
	DataManager *pDM, MBDynParser& HP);
extern MavlinkConsumer *ReadMarshConsumer(const std::string& name,
	DataManager *pDM, MBDynParser& HP);

#endif /* MODULE_MARSH_PUBSUB_H */
