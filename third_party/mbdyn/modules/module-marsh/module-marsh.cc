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

#include <cerrno>
#include <cstring>
#include <iostream>

#include <fcntl.h>
#include <netdb.h>
#include <sys/socket.h>
#include <unistd.h>

#include "dataman.h"
#include "userelem.h"

#include "module-marsh.h"

/* resolve host:port into a sockaddr_in; returns false on failure */
static bool
ResolveAddr(const char *host, unsigned short port, struct sockaddr_in& addr)
{
	struct addrinfo hints;
	memset(&hints, 0, sizeof(hints));
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_DGRAM;

	struct addrinfo *res = 0;
	if (getaddrinfo(host, 0, &hints, &res) != 0 || res == 0) {
		return false;
	}

	addr = *reinterpret_cast<struct sockaddr_in *>(res->ai_addr);
	addr.sin_port = htons(port);

	freeaddrinfo(res);

	return true;
}

ModuleMarsh::ModuleMarsh(
	unsigned uLabel_a, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: UserDefinedElem(uLabel_a, pDO),
m_pDM(pDM),
m_sock(-1),
m_systemId(1),
m_componentId(MARSH_COMP_ID_FLIGHT_MODEL),
m_heartbeatPeriod(1.),
m_dLastHeartbeatTime(-1.),
m_bManagerSeen(false),
m_dLastManagerSeenTime(-1.)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	marsh							\n"
"									\n"
"	Makes MBDyn act as a MARSH (https://marsh-sim.github.io/)	\n"
"	Flight Dynamics Module node: exchanges MAVLink messages	\n"
"	over UDP with the MARSH Manager.				\n"
"									\n"
"Syntax:								\n"
"	user defined: <label>, marsh,					\n"
"		[ manager address, (string) <addr>, ]			\n"
"		[ manager port, (int) <port>, ]			\n"
"		[ local port, (int) <port>, ]				\n"
"		[ system id, (int) <id>, ]				\n"
"		[ component id, (int) <id>, ]				\n"
"		[ heartbeat rate, (real) <Hz>, ]			\n"
"		{ send, (string) <message name>, ... ,			\n"
"		| subscribe, (string) <message name>, ... , }*		\n"
"	;								\n"
"									\n"
"	v1 messages: \"SIM_STATE\", \"MOTION_CUE_EXTRA\" (send),	\n"
"	\"MANUAL_CONTROL\" (subscribe); HEARTBEAT is automatic.	\n"
			<< std::endl);

		if (!HP.IsArg()) {
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	std::string managerAddr("127.0.0.1");
	if (HP.IsKeyWord("manager" "address")) {
		managerAddr = HP.GetStringWithDelims();
	}

	unsigned short managerPort = 24400;
	if (HP.IsKeyWord("manager" "port")) {
		managerPort = (unsigned short)HP.GetInt();
	}

	unsigned short localPort = 0;
	if (HP.IsKeyWord("local" "port")) {
		localPort = (unsigned short)HP.GetInt();
	}

	if (HP.IsKeyWord("system" "id")) {
		m_systemId = (uint8_t)HP.GetInt();
	}

	if (HP.IsKeyWord("component" "id")) {
		m_componentId = (uint8_t)HP.GetInt();
	}

	if (HP.IsKeyWord("heartbeat" "rate")) {
		doublereal dRate = HP.GetReal();
		if (dRate <= 0.) {
			silent_cerr("module-marsh: invalid heartbeat rate "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_heartbeatPeriod = 1./dRate;
	}

	for (;;) {
		if (HP.IsKeyWord("send")) {
			const std::string name = HP.GetStringWithDelims();
			MavlinkProducer *p = ReadMarshProducer(name, pDM, HP);
			if (!p) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			m_producers.push_back(p);

		} else if (HP.IsKeyWord("subscribe")) {
			const std::string name = HP.GetStringWithDelims();
			MavlinkConsumer *c = ReadMarshConsumer(name, pDM, HP);
			if (!c) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			m_consumers.push_back(c);

		} else {
			break;
		}
	}

	if (!ResolveAddr(managerAddr.c_str(), managerPort, m_managerAddr)) {
		silent_cerr("module-marsh: unable to resolve manager address \""
			<< managerAddr << "\" at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_sock = socket(AF_INET, SOCK_DGRAM, 0);
	if (m_sock < 0) {
		silent_cerr("module-marsh: socket() failed: "
			<< strerror(errno) << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	struct sockaddr_in localAddr;
	memset(&localAddr, 0, sizeof(localAddr));
	localAddr.sin_family = AF_INET;
	localAddr.sin_addr.s_addr = htonl(INADDR_ANY);
	localAddr.sin_port = htons(localPort);

	if (bind(m_sock, reinterpret_cast<struct sockaddr *>(&localAddr),
		sizeof(localAddr)) < 0)
	{
		silent_cerr("module-marsh: bind() failed: "
			<< strerror(errno) << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	int flags = fcntl(m_sock, F_GETFL, 0);
	fcntl(m_sock, F_SETFL, flags | O_NONBLOCK);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

ModuleMarsh::~ModuleMarsh(void)
{
	if (m_sock >= 0) {
		close(m_sock);
	}

	for (std::vector<MavlinkProducer *>::iterator i = m_producers.begin();
		i != m_producers.end(); ++i)
	{
		delete *i;
	}

	for (std::vector<MavlinkConsumer *>::iterator i = m_consumers.begin();
		i != m_consumers.end(); ++i)
	{
		delete *i;
	}
}

void
ModuleMarsh::SendMessage(const mavlink_message_t& msg)
{
	uint8_t buf[MAVLINK_MAX_PACKET_LEN];
	uint16_t len = mavlink_msg_to_send_buffer(buf, &msg);

	sendto(m_sock, buf, len, 0,
		reinterpret_cast<const struct sockaddr *>(&m_managerAddr),
		sizeof(m_managerAddr));
}

void
ModuleMarsh::SendHeartbeat(doublereal dTime)
{
	if (m_dLastHeartbeatTime >= 0.
		&& dTime - m_dLastHeartbeatTime < m_heartbeatPeriod)
	{
		return;
	}

	m_dLastHeartbeatTime = dTime;

	mavlink_message_t msg;
	mavlink_msg_heartbeat_pack(m_systemId, m_componentId, &msg,
		MAV_TYPE_GENERIC, MAV_AUTOPILOT_INVALID, 0, 0, MAV_STATE_ACTIVE);

	SendMessage(msg);
}

void
ModuleMarsh::SendProducers(doublereal dTime)
{
	uint32_t time_boot_ms = (uint32_t)(dTime*1000.);

	for (std::vector<MavlinkProducer *>::const_iterator i = m_producers.begin();
		i != m_producers.end(); ++i)
	{
		Vec3 gravity(::Zero3);
		const StructNode *pNode = (*i)->pGetNode();
		if (pNode != 0) {
			bGetGravity(pNode->GetXCurr(), gravity);
		}

		mavlink_message_t msg;
		(*i)->Pack(msg, time_boot_ms, m_systemId, m_componentId, gravity);
		SendMessage(msg);
	}
}

void
ModuleMarsh::Dispatch(const mavlink_message_t& msg, doublereal dTime)
{
	if (msg.msgid == MAVLINK_MSG_ID_HEARTBEAT) {
		m_bManagerSeen = true;
		m_dLastManagerSeenTime = dTime;
		return;
	}

	for (std::vector<MavlinkConsumer *>::const_iterator i = m_consumers.begin();
		i != m_consumers.end(); ++i)
	{
		if ((*i)->GetMsgId() == msg.msgid) {
			(*i)->Decode(msg);
		}
	}
}

void
ModuleMarsh::ReceivePending(doublereal dTime)
{
	uint8_t buf[512];
	for (;;) {
		ssize_t n = recvfrom(m_sock, buf, sizeof(buf), 0, 0, 0);
		if (n < 0) {
			/* EAGAIN/EWOULDBLOCK: nothing pending */
			break;
		}

		mavlink_message_t msg;
		mavlink_status_t status;
		for (ssize_t j = 0; j < n; j++) {
			if (mavlink_parse_char(MAVLINK_COMM_0, buf[j], &msg, &status)) {
				Dispatch(msg, dTime);
			}
		}
	}
}

void
ModuleMarsh::AfterConvergence(const VectorHandler& /* X */,
	const VectorHandler& /* XP */)
{
	doublereal dTime = m_pDM->dGetTime();

	ReceivePending(dTime);
	SendHeartbeat(dTime);
	SendProducers(dTime);
}

void
ModuleMarsh::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << GetLabel()
			<< " " << (m_bManagerSeen ? 1 : 0)
			<< std::endl;
	}
}

void
ModuleMarsh::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleMarsh::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal /* dCoef */,
	const VectorHandler& /* XCurr */,
	const VectorHandler& /* XPrimeCurr */)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
ModuleMarsh::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */,
	const VectorHandler& /* XCurr */,
	const VectorHandler& /* XPrimeCurr */)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}

bool
ModuleMarsh::GetPrivDataOwner(unsigned int i, MavlinkConsumer *& pConsumer,
	unsigned int& iLocal) const
{
	unsigned int offset = 0;

	for (std::vector<MavlinkConsumer *>::const_iterator ci = m_consumers.begin();
		ci != m_consumers.end(); ++ci)
	{
		unsigned int n = (*ci)->iGetNumPrivData();
		if (i > offset && i <= offset + n) {
			pConsumer = *ci;
			iLocal = i - offset;
			return true;
		}

		offset += n;
	}

	return false;
}

unsigned int
ModuleMarsh::iGetNumPrivData(void) const
{
	unsigned int n = 0;

	for (std::vector<MavlinkConsumer *>::const_iterator i = m_consumers.begin();
		i != m_consumers.end(); ++i)
	{
		n += (*i)->iGetNumPrivData();
	}

	return n;
}

unsigned int
ModuleMarsh::iGetPrivDataIdx(const char *s) const
{
	unsigned int offset = 0;

	for (std::vector<MavlinkConsumer *>::const_iterator i = m_consumers.begin();
		i != m_consumers.end(); ++i)
	{
		unsigned int idx = (*i)->iGetPrivDataIdx(s);
		if (idx > 0) {
			return offset + idx;
		}

		offset += (*i)->iGetNumPrivData();
	}

	return 0;
}

doublereal
ModuleMarsh::dGetPrivData(unsigned int i) const
{
	MavlinkConsumer *pConsumer;
	unsigned int iLocal;

	if (!GetPrivDataOwner(i, pConsumer, iLocal)) {
		silent_cerr("module-marsh(" << GetLabel() << "): "
			"private data " << i << " out of range" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pConsumer->dGetPrivData(iLocal);
}

int
ModuleMarsh::iGetNumConnectedNodes(void) const
{
	std::vector<const Node *> nodes;
	GetConnectedNodes(nodes);

	return nodes.size();
}

void
ModuleMarsh::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.clear();

	for (std::vector<MavlinkProducer *>::const_iterator i = m_producers.begin();
		i != m_producers.end(); ++i)
	{
		const StructNode *pNode = (*i)->pGetNode();
		if (pNode == 0) {
			continue;
		}

		bool bFound = false;
		for (std::vector<const Node *>::const_iterator j = connectedNodes.begin();
			j != connectedNodes.end(); ++j)
		{
			if (*j == pNode) {
				bFound = true;
				break;
			}
		}

		if (!bFound) {
			connectedNodes.push_back(pNode);
		}
	}
}

void
ModuleMarsh::SetValue(DataManager * /* pDM */,
	VectorHandler& /* X */, VectorHandler& /* XP */,
	SimulationEntity::Hints * /* ph */)
{
	NO_OP;
}

std::ostream&
ModuleMarsh::Restart(std::ostream& out) const
{
	return out << "# ModuleMarsh: not implemented" << std::endl;
}

unsigned int
ModuleMarsh::iGetInitialNumDof(void) const
{
	return 0;
}

void
ModuleMarsh::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleMarsh::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */)
{
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
ModuleMarsh::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */)
{
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<ModuleMarsh>;

	if (!SetUDE("marsh", rf)) {
		delete rf;

		silent_cerr("module-marsh: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}
