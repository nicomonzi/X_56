/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/rtai/pendulum-control/pendulum_control_rtai/pendulum_control.c,v 1.3 2007/02/26 00:34:50 masarati Exp $ */
/*
 * Real-Time Workshop code generation for Simulink model "pendulum_control.mdl".
 *
 * Model Version                        : 1.83
 * Real-Time Workshop file version      : 5.0 $Date: 2007/02/26 00:34:50 $
 * Real-Time Workshop file generated on : Thu Nov 27 14:52:23 2003
 * TLC version                          : 5.0 (Jun 18 2002)
 * C source code generated on           : Thu Nov 27 14:52:23 2003
 */

#include <math.h>
#include <string.h>
#include "pendulum_control.h"
#include "pendulum_control_private.h"

#include "mdl_info.h"

#include "pendulum_control_bio.c"

#include "pendulum_control_pt.c"

/* Block signals (auto storage) */
BlockIO rtB;

/* Block states (auto storage) */
D_Work rtDWork;

/* Parent Simstruct */
static rtModel_pendulum_control model_S;
rtModel_pendulum_control *const rtM_pendulum_control = &model_S;

/* Start for root system: '<Root>' */
void MdlStart(void)
{

  /* user code (Start function Header) */
  pendulum_control_InitializeParametersMap();

  /* Level2 S-Function Block: <Root>/S-Function1 (sfun_mbdyn_com_write) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[0];
    sfcnStart(rts);
  }

  /* Level2 S-Function Block: <Root>/S-Function2 (sfun_mbdyn_com_read) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[1];
    sfcnStart(rts);
  }

  /* Level2 S-Function Block: <Root>/C (sfun_rtai_scope) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[2];
    sfcnStart(rts);
  }

  /* Level2 S-Function Block: <Root>/err (sfun_rtai_scope) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[3];
    sfcnStart(rts);
  }

  /* Level2 S-Function Block: <Root>/omega (sfun_rtai_scope) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[4];
    sfcnStart(rts);
  }

  /* Level2 S-Function Block: <Root>/theta (sfun_rtai_scope) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[5];
    sfcnStart(rts);
  }
}

/* Outputs for root system: '<Root>' */
void MdlOutputs(int_T tid)
{
  /* tid is required for a uniform function interface. This system
   * is single rate, and in this case, tid is not accessed. */
  UNUSED_PARAMETER(tid);

  /* Level2 S-Function Block: <Root>/S-Function1 (sfun_mbdyn_com_write) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[0];
    sfcnOutputs(rts, tid);
  }

  /* Level2 S-Function Block: <Root>/S-Function2 (sfun_mbdyn_com_read) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[1];
    sfcnOutputs(rts, tid);
  }

  /* Sum: '<Root>/Sum2' incorporates:
   *   Constant: '<Root>/Constant1'
   */
  rtB.Sum2 = rtB.S_Function2_o1 - rtP.Constant1_Value;

  /* Sum: '<Root>/Sum' incorporates:
   *   Gain: '<Root>/Gain'
   *   Gain: '<Root>/Gain1'
   *
   * Regarding '<Root>/Gain':
   *   Gain value: rtP.Gain_Gain
   *
   * Regarding '<Root>/Gain1':
   *   Gain value: rtP.Gain1_Gain
   */
  rtB.Sum = - (rtB.Sum2 * rtP.Gain_Gain)
    - (rtB.S_Function2_o2 * rtP.Gain1_Gain);

  /* Level2 S-Function Block: <Root>/C (sfun_rtai_scope) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[2];
    sfcnOutputs(rts, tid);
  }

  /* Level2 S-Function Block: <Root>/err (sfun_rtai_scope) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[3];
    sfcnOutputs(rts, tid);
  }

  /* Level2 S-Function Block: <Root>/omega (sfun_rtai_scope) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[4];
    sfcnOutputs(rts, tid);
  }

  /* Level2 S-Function Block: <Root>/theta (sfun_rtai_scope) */
  {
    SimStruct *rts = rtM_pendulum_control->childSfunctions[5];
    sfcnOutputs(rts, tid);
  }

  /* Sum: '<Root>/Sum1' incorporates:
   *   Constant: '<Root>/Constant'
   */
  rtB.Sum1 = rtB.Sum - rtP.Constant_Value;
}

/* Update for root system: '<Root>' */
void MdlUpdate(int_T tid)
{

  /* tid is required for a uniform function interface. This system
   * is single rate, and in this case, tid is not accessed. */
  UNUSED_PARAMETER(tid);
}

/* Terminate for root system: '<Root>' */
void MdlTerminate(void)
{
  if(rtM_pendulum_control != NULL) {

    /* Level2 S-Function Block: <Root>/S-Function1 (sfun_mbdyn_com_write) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[0];
      sfcnTerminate(rts);
    }

    /* Level2 S-Function Block: <Root>/S-Function2 (sfun_mbdyn_com_read) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[1];
      sfcnTerminate(rts);
    }

    /* Level2 S-Function Block: <Root>/C (sfun_rtai_scope) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[2];
      sfcnTerminate(rts);
    }

    /* Level2 S-Function Block: <Root>/err (sfun_rtai_scope) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[3];
      sfcnTerminate(rts);
    }

    /* Level2 S-Function Block: <Root>/omega (sfun_rtai_scope) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[4];
      sfcnTerminate(rts);
    }

    /* Level2 S-Function Block: <Root>/theta (sfun_rtai_scope) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[5];
      sfcnTerminate(rts);
    }
  }
}

/* Helper function to make function calls from non-inlined S-functions */
int_T rt_CallSys(SimStruct *S, int_T element, int_T tid)
{
  (*(S)->callSys.fcns[element])((S)->callSys.args1[element],
   (S)->callSys.args2[element], tid);

  if (ssGetErrorStatus(S) != NULL) {
    return 0;
  } else {
    return 1;
  }
}

/* Function to initialize sizes */
void MdlInitializeSizes(void)
{
  rtM_pendulum_control->Sizes.numContStates = (0); /* Number of continuous states */
  rtM_pendulum_control->Sizes.numY = (0); /* Number of model outputs */
  rtM_pendulum_control->Sizes.numU = (0); /* Number of model inputs */
  rtM_pendulum_control->Sizes.sysDirFeedThru = (0); /* The model is not direct feedthrough */
  rtM_pendulum_control->Sizes.numSampTimes = (1); /* Number of sample times */
  rtM_pendulum_control->Sizes.numBlocks = (13); /* Number of blocks */
  rtM_pendulum_control->Sizes.numBlockIO = (6); /* Number of block outputs */
  rtM_pendulum_control->Sizes.numBlockPrms = (94); /* Sum of parameter "widths" */
}

/* Function to initialize sample times */
void MdlInitializeSampleTimes(void)
{
  /* task periods */
  rtM_pendulum_control->Timing.sampleTimes[0] = (0.001);

  /* task offsets */
  rtM_pendulum_control->Timing.offsetTimes[0] = (0.0);
}

/* Function to register the model */
rtModel_pendulum_control *pendulum_control(void)
{
  (void)memset((char *)rtM_pendulum_control, 0,
   sizeof(rtModel_pendulum_control));

  {
    /* Setup solver object */
    static RTWSolverInfo rt_SolverInfo;
    rtM_pendulum_control->solverInfo = (&rt_SolverInfo);

    rtsiSetSimTimeStepPtr(rtM_pendulum_control->solverInfo,
     &rtM_pendulum_control->Timing.simTimeStep);
    rtsiSetTPtr(rtM_pendulum_control->solverInfo,
     &rtmGetTPtr(rtM_pendulum_control));
    rtsiSetStepSizePtr(rtM_pendulum_control->solverInfo,
     &rtM_pendulum_control->Timing.stepSize);
    rtsiSetdXPtr(rtM_pendulum_control->solverInfo,
     &rtM_pendulum_control->ModelData.derivs);
    rtsiSetContStatesPtr(rtM_pendulum_control->solverInfo,
     &rtM_pendulum_control->ModelData.contStates);
    rtsiSetNumContStatesPtr(rtM_pendulum_control->solverInfo,
     &rtM_pendulum_control->Sizes.numContStates);
    rtsiSetErrorStatusPtr(rtM_pendulum_control->solverInfo,
     &rtmGetErrorStatus(rtM_pendulum_control));

    rtsiSetRTModelPtr(rtM_pendulum_control->solverInfo, rtM_pendulum_control);
  }

  /* timing info */
  {
    static time_T mdlPeriod[NSAMPLE_TIMES];
    static time_T mdlOffset[NSAMPLE_TIMES];
    static time_T mdlTaskTimes[NSAMPLE_TIMES];
    static int_T mdlTsMap[NSAMPLE_TIMES];
    static int_T mdlSampleHits[NSAMPLE_TIMES];

    {
      int_T i;

      for(i = 0; i < NSAMPLE_TIMES; i++) {
        mdlPeriod[i] = 0.0;
        mdlOffset[i] = 0.0;
        mdlTaskTimes[i] = 0.0;
      }
    }
    (void)memset((char_T *)&mdlTsMap[0], 0, 1 * sizeof(int_T));
    (void)memset((char_T *)&mdlSampleHits[0], 0, 1 * sizeof(int_T));

    rtM_pendulum_control->Timing.sampleTimes = (&mdlPeriod[0]);
    rtM_pendulum_control->Timing.offsetTimes = (&mdlOffset[0]);
    rtM_pendulum_control->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    rtmSetTPtr(rtM_pendulum_control, &mdlTaskTimes[0]);
    rtM_pendulum_control->Timing.sampleHits = (&mdlSampleHits[0]);
  }
  rtsiSetSolverMode(rtM_pendulum_control->solverInfo, SOLVER_MODE_SINGLETASKING);

  /*
   * initialize model vectors and cache them in SimStruct
   */

  /* block I/O */
  {
    void *b = (void *) &rtB;
    rtM_pendulum_control->ModelData.blockIO = (b);

    {
      int_T i;

      b =&rtB.S_Function1;
      for (i = 0; i < 6; i++) {
        ((real_T*)b)[i] = 0.0;
      }
    }
  }

  /* parameters */
  rtM_pendulum_control->ModelData.defaultParam = ((real_T *) &rtP);

  /* data type work */
  {
    void *dwork = (void *) &rtDWork;
    rtM_pendulum_control->Work.dwork = (dwork);
    (void)memset((char_T *) dwork, 0, sizeof(D_Work));
  }

  /* C API for Parameter Tuning and/or Signal Monitoring */
  {
    static ModelMappingInfo mapInfo;

    memset((char_T *) &mapInfo, 0, sizeof(mapInfo));

    /* block signal monitoring map */
    mapInfo.Signals.blockIOSignals = &rtBIOSignals[0];
    mapInfo.Signals.numBlockIOSignals = 6;

    /* parameter tuning maps */
    mapInfo.Parameters.blockTuning = &rtBlockTuning[0];
    mapInfo.Parameters.variableTuning = &rtVariableTuning[0];
    mapInfo.Parameters.parametersMap = rtParametersMap;
    mapInfo.Parameters.dimensionsMap = rtDimensionsMap;
    mapInfo.Parameters.numBlockTuning = 22;
    mapInfo.Parameters.numVariableTuning = 0;

    rtM_pendulum_control->SpecialInfo.mappingInfo = (&mapInfo);
  }

  /* Model specific registration */

  rtM_pendulum_control->modelName = ("pendulum_control");
  rtM_pendulum_control->path = ("pendulum_control");

  rtmSetTStart(rtM_pendulum_control, 0.0);
  rtM_pendulum_control->Timing.tFinal = (10.0);
  rtM_pendulum_control->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(rtM_pendulum_control->solverInfo, 0.001);
  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;

    rtM_pendulum_control->rtwLogInfo = (&rt_DataLoggingInfo);

    rtliSetLogFormat(rtM_pendulum_control->rtwLogInfo, 0);

    rtliSetLogMaxRows(rtM_pendulum_control->rtwLogInfo, 1000);

    rtliSetLogDecimation(rtM_pendulum_control->rtwLogInfo, 1);

    rtliSetLogVarNameModifier(rtM_pendulum_control->rtwLogInfo, "rt_");

    rtliSetLogT(rtM_pendulum_control->rtwLogInfo, "tout");

    rtliSetLogX(rtM_pendulum_control->rtwLogInfo, "");

    rtliSetLogXFinal(rtM_pendulum_control->rtwLogInfo, "");

    rtliSetLogXSignalInfo(rtM_pendulum_control->rtwLogInfo, NULL);

    rtliSetLogXSignalPtrs(rtM_pendulum_control->rtwLogInfo, NULL);

    rtliSetLogY(rtM_pendulum_control->rtwLogInfo, "");

    rtliSetLogYSignalInfo(rtM_pendulum_control->rtwLogInfo, NULL);

    rtliSetLogYSignalPtrs(rtM_pendulum_control->rtwLogInfo, NULL);
  }

  rtM_pendulum_control->Sizes.checksums[0] = (396058510U);
  rtM_pendulum_control->Sizes.checksums[1] = (790027248U);
  rtM_pendulum_control->Sizes.checksums[2] = (4068364517U);
  rtM_pendulum_control->Sizes.checksums[3] = (2378495770U);

  /* child S-Function registration */
  {
    static RTWSfcnInfo _sfcnInfo;
    RTWSfcnInfo *sfcnInfo = &_sfcnInfo;

    rtM_pendulum_control->sfcnInfo = (sfcnInfo);

    rtssSetErrorStatusPtr(sfcnInfo, &rtmGetErrorStatus(rtM_pendulum_control));
    rtssSetNumRootSampTimesPtr(sfcnInfo,
     &rtM_pendulum_control->Sizes.numSampTimes);
    rtssSetTPtrPtr(sfcnInfo, &rtmGetTPtr(rtM_pendulum_control));
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(rtM_pendulum_control));
    rtssSetTFinalPtr(sfcnInfo, &rtM_pendulum_control->Timing.tFinal);
    rtssSetTimeOfLastOutputPtr(sfcnInfo,
     &rtM_pendulum_control->Timing.timeOfLastOutput);
    rtssSetStepSizePtr(sfcnInfo, &rtM_pendulum_control->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo,
     &rtM_pendulum_control->Timing.stopRequestedFlag);
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
     &rtM_pendulum_control->ModelData.derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
     &rtM_pendulum_control->ModelData.zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo,
     &rtM_pendulum_control->ModelData.blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo, &rtM_pendulum_control->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
     &rtM_pendulum_control->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &rtM_pendulum_control->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &rtM_pendulum_control->solverInfo);
  }

  rtM_pendulum_control->Sizes.numSFcns = (6);

  /* register each child */
  {
    static SimStruct childSFunctions[6];
    static SimStruct *childSFunctionPtrs[6];

    (void)memset((char_T *)&childSFunctions[0], 0, sizeof(childSFunctions));
    rtM_pendulum_control->childSfunctions = (&childSFunctionPtrs[0]);
    {
      int_T i;

      for(i = 0; i < 6; i++) {
        rtM_pendulum_control->childSfunctions[i] = (&childSFunctions[i]);
      }
    }

    /* Level2 S-Function Block: pendulum_control/<Root>/S-Function1 (sfun_mbdyn_com_write) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[0];
      /* timing info */
      static time_T sfcnPeriod[1];
      static time_T sfcnOffset[1];
      static int_T sfcnTsMap[1];

      {
        int_T i;

        for(i = 0; i < 1; i++) {
          sfcnPeriod[i] = sfcnOffset[i] = 0.0;
        }
      }
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        static struct _ssBlkInfo2 _blkInfo2;
        struct _ssBlkInfo2 *blkInfo2 = &_blkInfo2;
        ssSetBlkInfo2Ptr(rts, blkInfo2);
        ssSetRTWSfcnInfo(rts, rtM_pendulum_control->sfcnInfo);
      }

      /* Allocate memory of model methods 2 */
      {
        static struct _ssSFcnModelMethods2 methods2;
        ssSetModelMethods2(rts, &methods2);
      }

      /* inputs */
      {
        static struct _ssPortInputs inputPortInfo[1];

        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts, &inputPortInfo[0]);

        /* port 0 */
        {

          static real_T const *sfcnUPtrs[1];
          sfcnUPtrs[0] = &rtB.Sum1;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }
      }

      /* outputs */
      {
        static struct _ssPortOutputs outputPortInfo[1];
        _ssSetNumOutputPorts(rts, 1);
        ssSetPortInfoForOutputs(rts, &outputPortInfo[0]);
        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &rtB.S_Function1));
        }
      }
      /* path info */
      ssSetModelName(rts, "S-Function1");
      ssSetPath(rts, "pendulum_control/S-Function1");
      ssSetRTModel(rts,rtM_pendulum_control);
      ssSetParentSS(rts, NULL);
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        static mxArray *sfcnParams[5];

        ssSetSFcnParamsCount(rts, 5);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);

        ssSetSFcnParam(rts, 0, &rtP.S_Function1_P1_Size[0]);
        ssSetSFcnParam(rts, 1, &rtP.S_Function1_P2_Size[0]);
        ssSetSFcnParam(rts, 2, &rtP.S_Function1_P3_Size[0]);
        ssSetSFcnParam(rts, 3, &rtP.S_Function1_P4_Size[0]);
        ssSetSFcnParam(rts, 4, &rtP.S_Function1_P5_Size[0]);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *) &rtDWork.S_Function1_IWORK[0]);
      ssSetPWork(rts, (void **) &rtDWork.S_Function1_PWORK[0]);
      {
        static struct _ssDWorkRecord dWorkRecord[2];

        ssSetSFcnDWork(rts, dWorkRecord);
        _ssSetNumDWork(rts, 2);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 2);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &rtDWork.S_Function1_IWORK[0]);

        /* PWORK */
        ssSetDWorkWidth(rts, 1, 2);
        ssSetDWorkDataType(rts, 1,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &rtDWork.S_Function1_PWORK[0]);
      }

      /* registration */
      sfun_mbdyn_com_write(rts);

      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.001);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */

      ssSetNumNonsampledZCs(rts, 0);
      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: pendulum_control/<Root>/S-Function2 (sfun_mbdyn_com_read) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[1];
      /* timing info */
      static time_T sfcnPeriod[1];
      static time_T sfcnOffset[1];
      static int_T sfcnTsMap[1];

      {
        int_T i;

        for(i = 0; i < 1; i++) {
          sfcnPeriod[i] = sfcnOffset[i] = 0.0;
        }
      }
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        static struct _ssBlkInfo2 _blkInfo2;
        struct _ssBlkInfo2 *blkInfo2 = &_blkInfo2;
        ssSetBlkInfo2Ptr(rts, blkInfo2);
        ssSetRTWSfcnInfo(rts, rtM_pendulum_control->sfcnInfo);
      }

      /* Allocate memory of model methods 2 */
      {
        static struct _ssSFcnModelMethods2 methods2;
        ssSetModelMethods2(rts, &methods2);
      }

      /* outputs */
      {
        static struct _ssPortOutputs outputPortInfo[2];
        _ssSetNumOutputPorts(rts, 2);
        ssSetPortInfoForOutputs(rts, &outputPortInfo[0]);
        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &rtB.S_Function2_o1));
        }
        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidth(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *) &rtB.S_Function2_o2));
        }
      }
      /* path info */
      ssSetModelName(rts, "S-Function2");
      ssSetPath(rts, "pendulum_control/S-Function2");
      ssSetRTModel(rts,rtM_pendulum_control);
      ssSetParentSS(rts, NULL);
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        static mxArray *sfcnParams[5];

        ssSetSFcnParamsCount(rts, 5);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);

        ssSetSFcnParam(rts, 0, &rtP.S_Function2_P1_Size[0]);
        ssSetSFcnParam(rts, 1, &rtP.S_Function2_P2_Size[0]);
        ssSetSFcnParam(rts, 2, &rtP.S_Function2_P3_Size[0]);
        ssSetSFcnParam(rts, 3, &rtP.S_Function2_P4_Size[0]);
        ssSetSFcnParam(rts, 4, &rtP.S_Function2_P5_Size[0]);
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *) &rtDWork.S_Function2_IWORK);
      ssSetPWork(rts, (void **) &rtDWork.S_Function2_PWORK);
      {
        static struct _ssDWorkRecord dWorkRecord[2];

        ssSetSFcnDWork(rts, dWorkRecord);
        _ssSetNumDWork(rts, 2);

        /* IWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &rtDWork.S_Function2_IWORK);

        /* PWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &rtDWork.S_Function2_PWORK);
      }

      /* registration */
      sfun_mbdyn_com_read(rts);

      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.001);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */

      ssSetNumNonsampledZCs(rts, 0);
      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);
      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: pendulum_control/<Root>/C (sfun_rtai_scope) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[2];
      /* timing info */
      static time_T sfcnPeriod[1];
      static time_T sfcnOffset[1];
      static int_T sfcnTsMap[1];

      {
        int_T i;

        for(i = 0; i < 1; i++) {
          sfcnPeriod[i] = sfcnOffset[i] = 0.0;
        }
      }
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        static struct _ssBlkInfo2 _blkInfo2;
        struct _ssBlkInfo2 *blkInfo2 = &_blkInfo2;
        ssSetBlkInfo2Ptr(rts, blkInfo2);
        ssSetRTWSfcnInfo(rts, rtM_pendulum_control->sfcnInfo);
      }

      /* Allocate memory of model methods 2 */
      {
        static struct _ssSFcnModelMethods2 methods2;
        ssSetModelMethods2(rts, &methods2);
      }

      /* inputs */
      {
        static struct _ssPortInputs inputPortInfo[2];

        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts, &inputPortInfo[0]);

        /* port 0 */
        {

          static real_T const *sfcnUPtrs[1];
          sfcnUPtrs[0] = &rtB.S_Function1;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {

          static real_T const *sfcnUPtrs[1];
          sfcnUPtrs[0] = &rtB.Sum;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* path info */
      ssSetModelName(rts, "C");
      ssSetPath(rts, "pendulum_control/C");
      ssSetRTModel(rts,rtM_pendulum_control);
      ssSetParentSS(rts, NULL);
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        static mxArray *sfcnParams[2];

        ssSetSFcnParamsCount(rts, 2);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);

        ssSetSFcnParam(rts, 0, &rtP.C_P1_Size[0]);
        ssSetSFcnParam(rts, 1, &rtP.C_P2_Size[0]);
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &rtDWork.C_PWORK);
      {
        static struct _ssDWorkRecord dWorkRecord[1];

        ssSetSFcnDWork(rts, dWorkRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &rtDWork.C_PWORK);
      }

      /* registration */
      sfun_rtai_scope(rts);

      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.001);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */

      ssSetNumNonsampledZCs(rts, 0);
      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: pendulum_control/<Root>/err (sfun_rtai_scope) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[3];
      /* timing info */
      static time_T sfcnPeriod[1];
      static time_T sfcnOffset[1];
      static int_T sfcnTsMap[1];

      {
        int_T i;

        for(i = 0; i < 1; i++) {
          sfcnPeriod[i] = sfcnOffset[i] = 0.0;
        }
      }
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        static struct _ssBlkInfo2 _blkInfo2;
        struct _ssBlkInfo2 *blkInfo2 = &_blkInfo2;
        ssSetBlkInfo2Ptr(rts, blkInfo2);
        ssSetRTWSfcnInfo(rts, rtM_pendulum_control->sfcnInfo);
      }

      /* Allocate memory of model methods 2 */
      {
        static struct _ssSFcnModelMethods2 methods2;
        ssSetModelMethods2(rts, &methods2);
      }

      /* inputs */
      {
        static struct _ssPortInputs inputPortInfo[1];

        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts, &inputPortInfo[0]);

        /* port 0 */
        {

          static real_T const *sfcnUPtrs[1];
          sfcnUPtrs[0] = &rtB.Sum2;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }
      }

      /* path info */
      ssSetModelName(rts, "err");
      ssSetPath(rts, "pendulum_control/err");
      ssSetRTModel(rts,rtM_pendulum_control);
      ssSetParentSS(rts, NULL);
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        static mxArray *sfcnParams[2];

        ssSetSFcnParamsCount(rts, 2);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);

        ssSetSFcnParam(rts, 0, &rtP.err_P1_Size[0]);
        ssSetSFcnParam(rts, 1, &rtP.err_P2_Size[0]);
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &rtDWork.err_PWORK);
      {
        static struct _ssDWorkRecord dWorkRecord[1];

        ssSetSFcnDWork(rts, dWorkRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &rtDWork.err_PWORK);
      }

      /* registration */
      sfun_rtai_scope(rts);

      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.001);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */

      ssSetNumNonsampledZCs(rts, 0);
      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: pendulum_control/<Root>/omega (sfun_rtai_scope) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[4];
      /* timing info */
      static time_T sfcnPeriod[1];
      static time_T sfcnOffset[1];
      static int_T sfcnTsMap[1];

      {
        int_T i;

        for(i = 0; i < 1; i++) {
          sfcnPeriod[i] = sfcnOffset[i] = 0.0;
        }
      }
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        static struct _ssBlkInfo2 _blkInfo2;
        struct _ssBlkInfo2 *blkInfo2 = &_blkInfo2;
        ssSetBlkInfo2Ptr(rts, blkInfo2);
        ssSetRTWSfcnInfo(rts, rtM_pendulum_control->sfcnInfo);
      }

      /* Allocate memory of model methods 2 */
      {
        static struct _ssSFcnModelMethods2 methods2;
        ssSetModelMethods2(rts, &methods2);
      }

      /* inputs */
      {
        static struct _ssPortInputs inputPortInfo[1];

        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts, &inputPortInfo[0]);

        /* port 0 */
        {

          static real_T const *sfcnUPtrs[1];
          sfcnUPtrs[0] = &rtB.S_Function2_o2;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }
      }

      /* path info */
      ssSetModelName(rts, "omega");
      ssSetPath(rts, "pendulum_control/omega");
      ssSetRTModel(rts,rtM_pendulum_control);
      ssSetParentSS(rts, NULL);
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        static mxArray *sfcnParams[2];

        ssSetSFcnParamsCount(rts, 2);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);

        ssSetSFcnParam(rts, 0, &rtP.omega_P1_Size[0]);
        ssSetSFcnParam(rts, 1, &rtP.omega_P2_Size[0]);
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &rtDWork.omega_PWORK);
      {
        static struct _ssDWorkRecord dWorkRecord[1];

        ssSetSFcnDWork(rts, dWorkRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &rtDWork.omega_PWORK);
      }

      /* registration */
      sfun_rtai_scope(rts);

      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.001);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */

      ssSetNumNonsampledZCs(rts, 0);
      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: pendulum_control/<Root>/theta (sfun_rtai_scope) */
    {
      SimStruct *rts = rtM_pendulum_control->childSfunctions[5];
      /* timing info */
      static time_T sfcnPeriod[1];
      static time_T sfcnOffset[1];
      static int_T sfcnTsMap[1];

      {
        int_T i;

        for(i = 0; i < 1; i++) {
          sfcnPeriod[i] = sfcnOffset[i] = 0.0;
        }
      }
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        static struct _ssBlkInfo2 _blkInfo2;
        struct _ssBlkInfo2 *blkInfo2 = &_blkInfo2;
        ssSetBlkInfo2Ptr(rts, blkInfo2);
        ssSetRTWSfcnInfo(rts, rtM_pendulum_control->sfcnInfo);
      }

      /* Allocate memory of model methods 2 */
      {
        static struct _ssSFcnModelMethods2 methods2;
        ssSetModelMethods2(rts, &methods2);
      }

      /* inputs */
      {
        static struct _ssPortInputs inputPortInfo[1];

        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts, &inputPortInfo[0]);

        /* port 0 */
        {

          static real_T const *sfcnUPtrs[1];
          sfcnUPtrs[0] = &rtB.S_Function2_o1;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }
      }

      /* path info */
      ssSetModelName(rts, "theta");
      ssSetPath(rts, "pendulum_control/theta");
      ssSetRTModel(rts,rtM_pendulum_control);
      ssSetParentSS(rts, NULL);
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        static mxArray *sfcnParams[2];

        ssSetSFcnParamsCount(rts, 2);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);

        ssSetSFcnParam(rts, 0, &rtP.theta_P1_Size[0]);
        ssSetSFcnParam(rts, 1, &rtP.theta_P2_Size[0]);
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &rtDWork.theta_PWORK);
      {
        static struct _ssDWorkRecord dWorkRecord[1];

        ssSetSFcnDWork(rts, dWorkRecord);
        _ssSetNumDWork(rts, 1);

        /* PWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &rtDWork.theta_PWORK);
      }

      /* registration */
      sfun_rtai_scope(rts);

      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.001);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */

      ssSetNumNonsampledZCs(rts, 0);
      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }
  }

  return rtM_pendulum_control;
}

