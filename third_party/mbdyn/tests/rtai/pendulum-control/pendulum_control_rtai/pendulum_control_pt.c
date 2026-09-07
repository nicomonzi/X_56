/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/rtai/pendulum-control/pendulum_control_rtai/pendulum_control_pt.c,v 1.3 2007/02/26 00:34:50 masarati Exp $ */
/*
 * pendulum_control_pt.c
 *
 * Real-Time Workshop code generation for Simulink model "pendulum_control.mdl".
 *
 * Model Version                        : 1.83
 * Real-Time Workshop file version      : 5.0 $Date: 2007/02/26 00:34:50 $
 * Real-Time Workshop file generated on : Thu Nov 27 14:52:23 2003
 * TLC version                          : 5.0 (Jun 18 2002)
 * C source code generated on           : Thu Nov 27 14:52:23 2003
 */

#ifndef _PT_INFO_pendulum_control_
#define _PT_INFO_pendulum_control_

#include "pt_info.h"

/* Tunable block parameters */

static const BlockTuning rtBlockTuning[] = {

  /* blockName, parameterName,
   * class, nRows, nCols, nDims, dimsOffset, source, dataType, numInstances,
   * mapOffset
   */

  /* S-Function */
  {"pendulum_control/S-Function1", "P1",
    {rt_MATRIX_COL_MAJOR, 1, 9, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 0}
  },
  /* S-Function */
  {"pendulum_control/S-Function1", "P2",
    {rt_MATRIX_COL_MAJOR, 1, 6, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 1}
  },
  /* S-Function */
  {"pendulum_control/S-Function1", "P3",
    {rt_MATRIX_COL_MAJOR, 1, 6, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 2}
  },
  /* S-Function */
  {"pendulum_control/S-Function1", "P4",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 3}
  },
  /* S-Function */
  {"pendulum_control/S-Function1", "P5",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 4}
  },
  /* S-Function */
  {"pendulum_control/S-Function2", "P1",
    {rt_MATRIX_COL_MAJOR, 1, 9, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 5}
  },
  /* S-Function */
  {"pendulum_control/S-Function2", "P2",
    {rt_MATRIX_COL_MAJOR, 1, 6, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 6}
  },
  /* S-Function */
  {"pendulum_control/S-Function2", "P3",
    {rt_MATRIX_COL_MAJOR, 1, 6, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 7}
  },
  /* S-Function */
  {"pendulum_control/S-Function2", "P4",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 8}
  },
  /* S-Function */
  {"pendulum_control/S-Function2", "P5",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 9}
  },
  /* Constant */
  {"pendulum_control/Constant1", "Value",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 10}
  },
  /* Gain */
  {"pendulum_control/Gain", "Gain",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 11}
  },
  /* Gain */
  {"pendulum_control/Gain1", "Gain",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 12}
  },
  /* S-Function */
  {"pendulum_control/C", "P1",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 13}
  },
  /* S-Function */
  {"pendulum_control/C", "P2",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 14}
  },
  /* S-Function */
  {"pendulum_control/err", "P1",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 15}
  },
  /* S-Function */
  {"pendulum_control/err", "P2",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 16}
  },
  /* S-Function */
  {"pendulum_control/omega", "P1",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 17}
  },
  /* S-Function */
  {"pendulum_control/omega", "P2",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 18}
  },
  /* S-Function */
  {"pendulum_control/theta", "P1",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 19}
  },
  /* S-Function */
  {"pendulum_control/theta", "P2",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 20}
  },
  /* Constant */
  {"pendulum_control/Constant", "Value",
    {rt_SCALAR, 1, 1, 2, -1, rt_SL_PARAM, SS_DOUBLE, 1, 21}
  },
  {NULL, NULL,
    {(ParamClass)0, 0, 0, 0, 0, (ParamSource)0, 0, 0, 0}
  }
};

/* Tunable variable parameters */

static const VariableTuning rtVariableTuning[] = {

  /* variableName,
   * class, nRows, nCols, nDims, dimsOffset, source, dataType, numInstances,
   * mapOffset
   */

  {NULL,
    {(ParamClass)0, 0, 0, 0, 0, (ParamSource)0, 0, 0, 0}
  }
};

static void * rtParametersMap[22];

void pendulum_control_InitializeParametersMap(void) {
  rtParametersMap[0] = &rtP.S_Function1_P1[0]; /* 0 */
  rtParametersMap[1] = &rtP.S_Function1_P2[0]; /* 1 */
  rtParametersMap[2] = &rtP.S_Function1_P3[0]; /* 2 */
  rtParametersMap[3] = &rtP.S_Function1_P4; /* 3 */
  rtParametersMap[4] = &rtP.S_Function1_P5; /* 4 */
  rtParametersMap[5] = &rtP.S_Function2_P1[0]; /* 5 */
  rtParametersMap[6] = &rtP.S_Function2_P2[0]; /* 6 */
  rtParametersMap[7] = &rtP.S_Function2_P3[0]; /* 7 */
  rtParametersMap[8] = &rtP.S_Function2_P4; /* 8 */
  rtParametersMap[9] = &rtP.S_Function2_P5; /* 9 */
  rtParametersMap[10] = &rtP.Constant1_Value; /* 10 */
  rtParametersMap[11] = &rtP.Gain_Gain; /* 11 */
  rtParametersMap[12] = &rtP.Gain1_Gain; /* 12 */
  rtParametersMap[13] = &rtP.C_P1;      /* 13 */
  rtParametersMap[14] = &rtP.C_P2;      /* 14 */
  rtParametersMap[15] = &rtP.err_P1;    /* 15 */
  rtParametersMap[16] = &rtP.err_P2;    /* 16 */
  rtParametersMap[17] = &rtP.omega_P1; /* 17 */
  rtParametersMap[18] = &rtP.omega_P2; /* 18 */
  rtParametersMap[19] = &rtP.theta_P1; /* 19 */
  rtParametersMap[20] = &rtP.theta_P2; /* 20 */
  rtParametersMap[21] = &rtP.Constant_Value; /* 21 */
}

static uint_T const rtDimensionsMap[] = {
  0                                     /* Dummy */
};

#endif                                  /* _PT_INFO_pendulum_control_ */
