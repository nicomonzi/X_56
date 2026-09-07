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

#ifndef MODULE_MARSH_H
#define MODULE_MARSH_H

#include <vector>
#include <netinet/in.h>

#include "mavlink_pubsub.h"

/*
 * MBDyn as a MARSH (https://marsh-sim.github.io/) Flight Dynamics Module
 * node: exchanges MAVLink messages over UDP with the MARSH Manager.
 *
 * Always sends/receives HEARTBEAT (session liveness); which other
 * messages are sent/subscribed to, and what MBDyn state feeds them, is
 * declared per-instance via `send`/`subscribe` clauses in the .mbd deck
 * (see mavlink_pubsub.h for the registry those clauses look up).
 */
class ModuleMarsh
: public UserDefinedElem {
private:
	DataManager *m_pDM;

	int m_sock;
	struct sockaddr_in m_managerAddr;

	uint8_t m_systemId;
	uint8_t m_componentId;

	doublereal m_heartbeatPeriod;
	doublereal m_dLastHeartbeatTime;

	bool m_bManagerSeen;
	doublereal m_dLastManagerSeenTime;

	std::vector<MavlinkProducer *> m_producers;
	std::vector<MavlinkConsumer *> m_consumers;

	void SendMessage(const mavlink_message_t& msg);
	void SendHeartbeat(doublereal dTime);
	void SendProducers(doublereal dTime);
	void ReceivePending(doublereal dTime);
	void Dispatch(const mavlink_message_t& msg, doublereal dTime);

	/* consumer owning global private-data index i (1-based), and its
	 * index local to that consumer */
	bool GetPrivDataOwner(unsigned int i, MavlinkConsumer *& pConsumer,
		unsigned int& iLocal) const;

public:
	ModuleMarsh(unsigned uLabel_a, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleMarsh(void);

	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	unsigned int iGetNumPrivData(void) const;
	unsigned int iGetPrivDataIdx(const char *s) const;
	doublereal dGetPrivData(unsigned int i) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	virtual void AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP);
};

#endif /* MODULE_MARSH_H */
