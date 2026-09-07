/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/rtai/pendulum-control/pendulum_control_rtai/pendulum_control_bio.c,v 1.3 2007/02/26 00:34:50 masarati Exp $ */
/*
 * pendulum_control_bio.c
 *
 * Real-Time Workshop code generation for Simulink model "pendulum_control.mdl".
 *
 * Model Version                        : 1.83
 * Real-Time Workshop file version      : 5.0 $Date: 2007/02/26 00:34:50 $
 * Real-Time Workshop file generated on : Thu Nov 27 14:52:23 2003
 * TLC version                          : 5.0 (Jun 18 2002)
 * C source code generated on           : Thu Nov 27 14:52:23 2003
 */

#ifndef BLOCK_IO_SIGNALS
#define BLOCK_IO_SIGNALS

#include "bio_sig.h"

/* Block output signal information */
static const BlockIOSignals rtBIOSignals[] = {

  /* blockName,
   * signalName, portNumber, signalWidth, signalAddr, dtName, dtSize
   */

  {
    "pendulum_control/S-Function1",
    "NULL", 0, 1, &rtB.S_Function1, "double", sizeof(real_T)
  },
  {
    "pendulum_control/S-Function2",
    "NULL", 0, 1, &rtB.S_Function2_o1, "double", sizeof(real_T)
  },
  {
    "pendulum_control/S-Function2",
    "NULL", 1, 1, &rtB.S_Function2_o2, "double", sizeof(real_T)
  },
  {
    "pendulum_control/Sum2",
    "NULL", 0, 1, &rtB.Sum2, "double", sizeof(real_T)
  },
  {
    "pendulum_control/Sum",
    "NULL", 0, 1, &rtB.Sum, "double", sizeof(real_T)
  },
  {
    "pendulum_control/Sum1",
    "NULL", 0, 1, &rtB.Sum1, "double", sizeof(real_T)
  },
  {
    NULL, NULL, 0, 0, 0, NULL, 0
  }
};

#endif                                  /* BLOCK_IO_SIGNALS */
