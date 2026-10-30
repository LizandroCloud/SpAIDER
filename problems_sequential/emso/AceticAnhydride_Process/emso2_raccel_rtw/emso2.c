#include <math.h>
#include "emso2.h"
#include "emso2_private.h"
#include "emso2_dt.h"

const int_T gblNumToFiles = 0;
const int_T gblNumFrFiles = 0;
const int_T gblNumFrWksBlocks = 0;

#ifndef MATLAB_ROOTQ
#  define MATLAB_ROOTQ                 QUOTE(MATLAB_ROOT)
#endif

char gbl_raccel_MATLABRoot[] = MATLAB_ROOTQ;

#ifdef RSIM_WITH_SOLVER_MULTITASKING

const boolean_T gbl_raccel_isMultitasking = 1;

#else

const boolean_T gbl_raccel_isMultitasking = 0;

#endif

const boolean_T gbl_raccel_tid01eq = 0;
const int_T gbl_raccel_NumST = 2;
const char_T *gbl_raccel_Version = "7.1  (R2008a)  23-Jan-2008";
void raccel_setup_MMIStateLog(SimStruct* S)
{

#ifdef UseMMIDataLogging

  rt_FillStateSigInfoFromMMI(ssGetRTWLogInfo(S),&ssGetErrorStatus(S));

#endif

}

const int_T gblNumRootInportBlks = 0;
const int_T gblNumModelInputs = 0;
extern rtInportTUtable *gblInportTUtables;
extern const char *gblInportFileName;
const int_T gblInportDataTypeIdx[] = { -1 };

const int_T gblInportDims[] = { -1 } ;

const int_T gblInportComplex[] = { -1 };

const int_T gblInportInterpoFlag[] = { -1 };

const int_T gblInportContinuous[] = { -1 };

#include "simstruc.h"
#include "fixedpoint.h"

BlockIO rtB;
ContinuousStates rtX;
D_Work rtDWork;
static SimStruct model_S;
SimStruct *const rtS = &model_S;
void MdlInitialize(void)
{
  {
    SimStruct *rts = ssGetSFunction(rtS, 0);
    sfcnInitializeConditions(rts);
    if (ssGetErrorStatus(rts) != NULL)
      return;
  }
}

void MdlStart(void)
{
  {
    SimStruct *rts = ssGetSFunction(rtS, 0);
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != NULL)
      return;
  }

  MdlInitialize();
}

void MdlOutputs(int_T tid)
{
  {
    SimStruct *rts = ssGetSFunction(rtS, 0);
    sfcnOutputs(rts, 0);
  }

  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  UNUSED_PARAMETER(tid);
}

void MdlDerivatives(void)
{
  {
    SimStruct *rts = ssGetSFunction(rtS, 0);
    real_T *sfcndX_fx = (real_T *) &((StateDerivatives *) ssGetdX(rtS))
      ->EMSOmodel_CSTATE[0];
    ssSetdX(rts, sfcndX_fx);
    sfcnDerivatives(rts);
    if (ssGetErrorStatus(rts) != NULL)
      return;
  }
}

void MdlProjection(void)
{
}

void MdlTerminate(void)
{
  {
    SimStruct *rts = ssGetSFunction(rtS, 0);
    sfcnTerminate(rts);
  }
}

void MdlInitializeSizes(void)
{
  ssSetNumContStates(rtS, 40);
  ssSetNumY(rtS, 0);
  ssSetNumU(rtS, 0);
  ssSetDirectFeedThrough(rtS, 0);
  ssSetNumSampleTimes(rtS, 2);
  ssSetNumBlocks(rtS, 3);
  ssSetNumBlockIO(rtS, 1);
  ssSetNumBlockParams(rtS, 19);
}

void MdlInitializeSampleTimes(void)
{
  ssSetSampleTime(rtS, 0, 0.0);
  ssSetSampleTime(rtS, 1, 0.0);
  ssSetOffsetTime(rtS, 0, 0.0);
  ssSetOffsetTime(rtS, 1, 1.0);
}

void raccel_set_checksum( SimStruct *rtS)
{
  ssSetChecksumVal(rtS, 0, 2422102931U);
  ssSetChecksumVal(rtS, 1, 1387238813U);
  ssSetChecksumVal(rtS, 2, 3087840110U);
  ssSetChecksumVal(rtS, 3, 3840014613U);
}

SimStruct *raccel_register_model(void)
{
  static struct _ssMdlInfo mdlInfo;
  (void) memset((char *)rtS,0,
                sizeof(SimStruct));
  (void) memset((char *)&mdlInfo,0,
                sizeof(struct _ssMdlInfo));
  ssSetMdlInfoPtr(rtS, &mdlInfo);

  {
    static time_T mdlPeriod[NSAMPLE_TIMES];
    static time_T mdlOffset[NSAMPLE_TIMES];
    static time_T mdlTaskTimes[NSAMPLE_TIMES];
    static int_T mdlTsMap[NSAMPLE_TIMES];
    static int_T mdlSampleHits[NSAMPLE_TIMES];
    static boolean_T mdlTNextWasAdjustedPtr[NSAMPLE_TIMES];
    static int_T mdlPerTaskSampleHits[NSAMPLE_TIMES * NSAMPLE_TIMES];

    {
      int_T i;
      for (i = 0; i < NSAMPLE_TIMES; i++) {
        mdlPeriod[i] = 0.0;
        mdlOffset[i] = 0.0;
        mdlTaskTimes[i] = 0.0;
        mdlTsMap[i] = i;
        mdlSampleHits[i] = 1;
      }
    }

    ssSetSampleTimePtr(rtS, &mdlPeriod[0]);
    ssSetOffsetTimePtr(rtS, &mdlOffset[0]);
    ssSetSampleTimeTaskIDPtr(rtS, &mdlTsMap[0]);
    ssSetTPtr(rtS, &mdlTaskTimes[0]);
    ssSetSampleHitPtr(rtS, &mdlSampleHits[0]);
    ssSetTNextWasAdjustedPtr(rtS, &mdlTNextWasAdjustedPtr[0]);
    ssSetPerTaskSampleHitsPtr(rtS, &mdlPerTaskSampleHits[0]);
  }

  ssSetSolverMode(rtS, SOLVER_MODE_SINGLETASKING);

  {
    ssSetBlockIO(rtS, ((void *) &rtB));

    {
      ((real_T*)&rtB.EMSOmodel)[0] = 0.0;
    }
  }

  ssSetDefaultParam(rtS, (real_T *) &rtP);

  {
    real_T *x = (real_T *) &rtX;
    ssSetContStates(rtS, x);
    (void) memset((void *)x,0,
                  sizeof(ContinuousStates));
  }

  {
    void *dwork = (void *) &rtDWork;
    ssSetRootDWork(rtS, dwork);
    (void) memset((void *) dwork,0,
                  sizeof(D_Work));
  }

  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo,0,
                  sizeof(dtInfo));
    ssSetModelMappingInfo(rtS, &dtInfo);
    dtInfo.numDataTypes = 14;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];
    dtInfo.B = &rtBTransTable;
    dtInfo.P = &rtPTransTable;
  }

  ssSetRootSS(rtS, rtS);
  ssSetVersion(rtS, SIMSTRUCT_VERSION_LEVEL2);
  ssSetModelName(rtS, "emso2");
  ssSetPath(rtS, "emso2");
  ssSetTStart(rtS, 0.0);
  ssSetTFinal(rtS, 10.0);

  {
    static RTWLogInfo rt_DataLoggingInfo;
    ssSetRTWLogInfo(rtS, &rt_DataLoggingInfo);
    rtliSetLogXSignalInfo(ssGetRTWLogInfo(rtS), NULL);
    rtliSetLogXSignalPtrs(ssGetRTWLogInfo(rtS), NULL);
    rtliSetLogT(ssGetRTWLogInfo(rtS), "tout");
    rtliSetLogX(ssGetRTWLogInfo(rtS), "");
    rtliSetLogXFinal(ssGetRTWLogInfo(rtS), "");
    rtliSetSigLog(ssGetRTWLogInfo(rtS), "");
    rtliSetLogVarNameModifier(ssGetRTWLogInfo(rtS), "none");
    rtliSetLogFormat(ssGetRTWLogInfo(rtS), 0);
    rtliSetLogMaxRows(ssGetRTWLogInfo(rtS), 1000);
    rtliSetLogDecimation(ssGetRTWLogInfo(rtS), 1);
    rtliSetLogY(ssGetRTWLogInfo(rtS), "");
    rtliSetLogYSignalInfo(ssGetRTWLogInfo(rtS), NULL);
    rtliSetLogYSignalPtrs(ssGetRTWLogInfo(rtS), NULL);
  }

  {
    static ssSolverInfo slvrInfo;
    static boolean_T contStatesDisabled[40];
    static real_T solverAbsTol[40] = { 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006,
      1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006,
      1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006,
      1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006,
      1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006,
      1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006, 1.0E-006,
      1.0E-006 };

    static boolean_T solverAutoAbsTol[40] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1 };

    ssSetSolverRelTol(rtS, 0.001);
    ssSetSolverAbsTol(rtS, solverAbsTol);
    ssSetSolverAutoAbsTol(rtS, solverAutoAbsTol);
    ssSetStepSize(rtS, 0.0);
    ssSetMinStepSize(rtS, 0.0);
    ssSetMaxNumMinSteps(rtS, -1);
    ssSetMinStepViolatedError(rtS, 0);
    ssSetMaxStepSize(rtS, 7.0);
    ssSetSolverMaxOrder(rtS, -1);
    ssSetSolverRefineFactor(rtS, 1);
    ssSetOutputTimes(rtS, NULL);
    ssSetNumOutputTimes(rtS, 0);
    ssSetOutputTimesOnly(rtS, 0);
    ssSetOutputTimesIndex(rtS, 0);
    ssSetZCCacheNeedsReset(rtS, 0);
    ssSetDerivCacheNeedsReset(rtS, 0);
    ssSetNumNonContDerivSigInfos(rtS, 0);
    ssSetNonContDerivSigInfos(rtS, NULL);
    ssSetSolverInfo(rtS, &slvrInfo);
    ssSetSolverName(rtS, "ode23s");
    ssSetVariableStepSolver(rtS, 1);
    ssSetSolverConsistencyChecking(rtS, 0);
    ssSetSolverAdaptiveZcDetection(rtS, 0);
    ssSetSolverRobustResetMethod(rtS, 0);
    ssSetSolverStateProjection(rtS, 0);
    ssSetSolverMassMatrixType(rtS, (ssMatrixType)0);
    ssSetSolverMassMatrixNzMax(rtS, 0);
    ssSetModelOutputs(rtS, MdlOutputs);
    ssSetModelLogData(rtS, rt_UpdateTXYLogVars);
    ssSetModelUpdate(rtS, MdlUpdate);
    ssSetModelDerivatives(rtS, MdlDerivatives);
    ssSetSolverMaxConsecutiveMinStep(rtS, 1);
    ssSetSolverShapePreserveControl(rtS, 2);
    ssSetTNextTid(rtS, INT_MIN);
    ssSetTNext(rtS, rtMinusInf);
    ssSetSolverNeedsReset(rtS);
    ssSetNumNonsampledZCs(rtS, 0);
    ssSetContStateDisabled(rtS, contStatesDisabled);
    ssSetSolverMaxConsecutiveMinStep(rtS, 1);
  }

  ssSetChecksumVal(rtS, 0, 2422102931U);
  ssSetChecksumVal(rtS, 1, 1387238813U);
  ssSetChecksumVal(rtS, 2, 3087840110U);
  ssSetChecksumVal(rtS, 3, 3840014613U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    ssSetRTWExtModeInfo(rtS, &rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(&rt_ExtModeInfo, &ssGetModelMappingInfo(rtS));
    rteiSetChecksumsPtr(&rt_ExtModeInfo, ssGetChecksums(rtS));
    rteiSetTPtr(&rt_ExtModeInfo, ssGetTPtr(rtS));
  }

  ssSetNumSFunctions(rtS, 1);

  {
    static SimStruct childSFunctions[1];
    static SimStruct *childSFunctionPtrs[1];
    (void) memset((void *)&childSFunctions[0],0,
                  sizeof(childSFunctions));
    ssSetSFunctions(rtS, &childSFunctionPtrs[0]);
    ssSetSFunction(rtS, 0, &childSFunctions[0]);

    {
      SimStruct *rts = ssGetSFunction(rtS, 0);
      static time_T sfcnPeriod[1];
      static time_T sfcnOffset[1];
      static int_T sfcnTsMap[1];
      (void) memset((void*)sfcnPeriod,0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset,0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);
      ssSetMdlInfoPtr(rts, ssGetMdlInfoPtr(rtS));

      {
        static struct _ssSFcnModelMethods2 methods2;
        ssSetModelMethods2(rts, &methods2);
      }

      {
        static struct _ssPortInputs inputPortInfo[1];
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts, &inputPortInfo[0]);

        {
          static real_T const *sfcnUPtrs[1];
          sfcnUPtrs[0] = &rtP.Constant1_Value;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }
      }

      {
        static struct _ssPortOutputs outputPortInfo[1];
        ssSetPortInfoForOutputs(rts, &outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &rtB.EMSOmodel));
        }
      }

      ssSetContStates(rts, &rtX.EMSOmodel_CSTATE[0]);
      ssSetModelName(rts, "EMSO model");
      ssSetPath(rts, "emso2/EMSO model");
      if (ssGetRTModel(rtS) == NULL) {
        ssSetParentSS(rts, rtS);
        ssSetRootSS(rts, ssGetRootSS(rtS));
      } else {
        ssSetRTModel(rts,ssGetRTModel(rtS));
        ssSetParentSS(rts, NULL);
        ssSetRootSS(rts, rts);
      }

      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      {
        static mxArray *sfcnParams[3];
        ssSetSFcnParamsCount(rts, 3);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)&rtP.EMSOmodel_P1_Size[0]);
        ssSetSFcnParam(rts, 1, (mxArray*)&rtP.EMSOmodel_P2_Size[0]);
        ssSetSFcnParam(rts, 2, (mxArray*)&rtP.EMSOmodel_P3_Size[0]);
      }

      ssSetIWork(rts, (int_T *) &rtDWork.EMSOmodel_IWORK[0]);
      ssSetPWork(rts, (void **) &rtDWork.EMSOmodel_PWORK[0]);

      {
        static struct _ssDWorkRecord dWorkRecord[2];
        static struct _ssDWorkAuxRecord dWorkAuxRecord[2];
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 2);
        ssSetDWorkWidth(rts, 0, 4);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &rtDWork.EMSOmodel_IWORK[0]);
        ssSetDWorkWidth(rts, 1, 13);
        ssSetDWorkDataType(rts, 1,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &rtDWork.EMSOmodel_PWORK[0]);
      }

      emso_sf_c(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);
      ssSetSampleTime(rts, 0, 0.0);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;
      ssSetNumNonsampledZCs(rts, 0);
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }
  }

  return rtS;
}
