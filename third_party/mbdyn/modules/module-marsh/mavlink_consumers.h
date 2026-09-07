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

#ifndef MODULE_MARSH_CONSUMERS_H
#define MODULE_MARSH_CONSUMERS_H

#include "mavlink_pubsub.h"

/*
 * Decodes MANUAL_CONTROL and exposes its axes/buttons as private data
 * "x", "y", "z", "r", "buttons" (consumable, once bound under a `subscribe,
 * "MANUAL_CONTROL"` clause, as e.g.
 * `element, <marsh label>, loadable, string, "MANUAL_CONTROL.x"`).
 */
class ManualControlConsumer : public MavlinkConsumer {
private:
	int16_t m_x, m_y, m_z, m_r;
	uint16_t m_buttons;

public:
	ManualControlConsumer(DataManager *pDM, MBDynParser& HP);

	uint32_t GetMsgId(void) const override;
	void Decode(const mavlink_message_t& msg) override;

	unsigned int iGetNumPrivData(void) const override;
	unsigned int iGetPrivDataIdx(const char *s) const override;
	doublereal dGetPrivData(unsigned int i) const override;
};

#endif /* MODULE_MARSH_CONSUMERS_H */
