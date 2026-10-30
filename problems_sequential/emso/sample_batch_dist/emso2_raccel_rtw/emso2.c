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
const int_T gbl_raccel_NumST = 1;
const char_T *gbl_raccel_Version = "7.1  (R2008a)  23-Jan-2008";
void raccel_setup_MMIStateLog(SimStruct* S)
{

#ifdef UseMMIDataLogging

  rt_FillStateSigInfoFromMMI(ssGetRTWLogInfo(S),&ssGetErrorStatus(S));

#endif

}

const int_T gblNumRootInportBlks = 1;
const int_T gblNumModelInputs = 1;
extern rtInportTUtable *gblInportTUtables;
extern const char *gblInportFileName;
const int_T gblInportDataTypeIdx[] = { 0 };

const int_T gblInportDims[] = { 1, 1 } ;

const int_T gblInportComplex[] = { 0 };

const int_T gblInportInterpoFlag[] = { 1 };

const int_T gblInportContinuous[] = { 0 };

#include "simstruc.h"
#include "fixedpoint.h"

BlockIO rtB;
D_Work rtDWork;
ExternalInputs rtU;
ExternalOutputs rtY;
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
  if (gblInportFileName != NULL) {
    real_T time = ssGetT(rtS);
    int_T currTimeIdx;
    int_T i;

    {
      int k = 1;
      if (gblInportTUtables[0].nTimePoints == 1) {
        k = 0;
      }

      currTimeIdx = rt_getTimeIdx(gblInportTUtables[0].time, time,
        gblInportTUtables[0].nTimePoints,
        gblInportTUtables[0].currTimeIdx,
        1);
      gblInportTUtables[0].currTimeIdx = currTimeIdx;
      for (i = 0; i < 1; i++) {
        real_T* realPtr1 = (real_T*)gblInportTUtables[0].ur +
          i*gblInportTUtables[0].nTimePoints +currTimeIdx;
        real_T* realPtr2 = realPtr1 + 1*k;
        (void)rt_Interpolate_Datatype(
          realPtr1,
          realPtr2,
          &rtU.In1,
          time,
          gblInportTUtables[0].time[currTimeIdx],
          gblInportTUtables[0].time[currTimeIdx + k],
          gblInportTUtables[0].uDataType );
      }
    }
  }

  if (ssIsSampleHit(rtS, 0, 0)) {
    {
      SimStruct *rts = ssGetSFunction(rtS, 0);
      sfcnOutputs(rts, 0);
    }

    rtY.Out1[0] = rtB.EMSOmodel[0];
    rtY.Out1[1] = rtB.EMSOmodel[1];
    rtY.Out1[2] = rtB.EMSOmodel[2];
    rtY.Out1[3] = rtB.EMSOmodel[3];
  }

  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  if (ssIsSampleHit(rtS, 0, 0)) {
  }

  UNUSED_PARAMETER(tid);
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
  ssSetNumContStates(rtS, 0);
  ssSetNumY(rtS, 4);
  ssSetNumU(rtS, 1);
  ssSetDirectFeedThrough(rtS, 1);
  ssSetNumSampleTimes(rtS, 1);
  ssSetNumBlocks(rtS, 8);
  ssSetNumBlockIO(rtS, 1);
  ssSetNumBlockParams(rtS, 26);
}

void MdlInitializeSampleTimes(void)
{
  ssSetSampleTime(rtS, 0, 0.1);
  ssSetOffsetTime(rtS, 0, 0.0);
}

void raccel_set_checksum( SimStruct *rtS)
{
  ssSetChecksumVal(rtS, 0, 810755476U);
  ssSetChecksumVal(rtS, 1, 2513967157U);
  ssSetChecksumVal(rtS, 2, 2290698586U);
  ssSetChecksumVal(rtS, 3, 2399357202U);
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
      int_T i;
      void *pVoidBlockIORegion;
      pVoidBlockIORegion = (void *)(&rtB.EMSOmodel[0]);
      for (i = 0; i < 4; i++) {
        ((real_T*)pVoidBlockIORegion)[i] = 0.0;
      }
    }
  }

  {
    ssSetU(rtS, ((void*) &rtU));
    rtU.In1 = 0.0;
  }

  {
    ssSetY(rtS, &rtY);

    {
      int_T i;
      for (i = 0; i < 4; i++) {
        rtY.Out1[i] = 0.0;
      }
    }
  }

  ssSetDefaultParam(rtS, (real_T *) &rtP);

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
  ssSetTFinal(rtS, 36.0);

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

    {
      static void * rt_LoggedOutputSignalPtrs[] = {
        &rtY.Out1[0]
      };

      rtliSetLogYSignalPtrs(ssGetRTWLogInfo(rtS), ((LogSignalPtrsType)
        rt_LoggedOutputSignalPtrs));
    }

    {
      static int_T rt_LoggedOutputWidths[] = {
        4
      };

      static int_T rt_LoggedOutputNumDimensions[] = {
        1
      };

      static int_T rt_LoggedOutputDimensions[] = {
        4
      };

      static boolean_T rt_LoggedOutputIsVarDims[] = {
        0
      };

      static int_T* rt_LoggedCurrentSignalDimensions[] = {
        NULL
      };

      static BuiltInDTypeId rt_LoggedOutputDataTypeIds[] = {
        SS_DOUBLE
      };

      static int_T rt_LoggedOutputComplexSignals[] = {
        0
      };

      static const char_T *rt_LoggedOutputLabels[] = {
        "" };

      static const char_T *rt_LoggedOutputBlockNames[] = {
        "emso2/Out1" };

      static RTWLogDataTypeConvert rt_RTWLogDataTypeConvert[] = {
        { 0, SS_DOUBLE, SS_DOUBLE, 0, 0, 0, 1.0, 0, 0.0 }
      };

      static RTWLogSignalInfo rt_LoggedOutputSignalInfo[] = {
        {
          1,
          rt_LoggedOutputWidths,
          rt_LoggedOutputNumDimensions,
          rt_LoggedOutputDimensions,
          rt_LoggedOutputIsVarDims,
          rt_LoggedCurrentSignalDimensions,
          rt_LoggedOutputDataTypeIds,
          rt_LoggedOutputComplexSignals,
          NULL,

          { rt_LoggedOutputLabels },
          NULL,
          NULL,
          NULL,

          { rt_LoggedOutputBlockNames },

          { NULL },
          NULL,
          rt_RTWLogDataTypeConvert
        }
      };

      rtliSetLogYSignalInfo(ssGetRTWLogInfo(rtS), rt_LoggedOutputSignalInfo);
      rt_LoggedCurrentSignalDimensions[0] = &rt_LoggedOutputWidths[0];
    }

    rtliSetLogY(ssGetRTWLogInfo(rtS), "yout");
  }

  {
    static ssSolverInfo slvrInfo;
    ssSetStepSize(rtS, 0.1);
    ssSetMinStepSize(rtS, 2.2250738585072014E-308);
    ssSetMaxNumMinSteps(rtS, -1);
    ssSetMinStepViolatedError(rtS, 0);
    ssSetMaxStepSize(rtS, 0.1);
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
    ssSetSolverName(rtS, "VariableStepDiscrete");
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
    ssSetTNextTid(rtS, INT_MIN);
    ssSetTNext(rtS, rtMinusInf);
    ssSetSolverNeedsReset(rtS);
    ssSetNumNonsampledZCs(rtS, 0);
  }

  ssSetChecksumVal(rtS, 0, 810755476U);
  ssSetChecksumVal(rtS, 1, 2513967157U);
  ssSetChecksumVal(rtS, 2, 2290698586U);
  ssSetChecksumVal(rtS, 3, 2399357202U);

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
          sfcnUPtrs[0] = &rtP.Constant3_Value;
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
          ssSetOutputPortWidth(rts, 0, 4);
          ssSetOutputPortSignal(rts, 0, ((real_T *) rtB.EMSOmodel));
        }
      }

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
        static mxArray *sfcnParams[4];
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)&rtP.EMSOmodel_P1_Size[0]);
        ssSetSFcnParam(rts, 1, (mxArray*)&rtP.EMSOmodel_P2_Size[0]);
        ssSetSFcnParam(rts, 2, (mxArray*)&rtP.EMSOmodel_P3_Size[0]);
        ssSetSFcnParam(rts, 3, (mxArray*)&rtP.EMSOmodel_P4_Size[0]);
      }

      ssSetIWork(rts, (int_T *) &rtDWork.EMSOmodel_IWORK[0]);
      ssSetPWork(rts, (void **) &rtDWork.EMSOmodel_PWORK[0]);

      {
        static struct _ssDWorkRecord dWorkRecord[2];
        static struct _ssDWorkAuxRecord dWorkAuxRecord[2];
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 2);
        ssSetDWorkWidth(rts, 0, 3);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &rtDWork.EMSOmodel_IWORK[0]);
        ssSetDWorkWidth(rts, 1, 13);
        ssSetDWorkDataType(rts, 1,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &rtDWork.EMSOmodel_PWORK[0]);
      }

      emso_sf_d(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);
      ssSetSampleTime(rts, 0, 0.1);
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
