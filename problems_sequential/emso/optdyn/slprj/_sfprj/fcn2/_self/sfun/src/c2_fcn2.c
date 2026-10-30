/* Include files */

#include "fcn2_sfun.h"
#include "c2_fcn2.h"

/* Type Definitions */

/* Named Constants */
#define c2_IN_NO_ACTIVE_CHILD          (0)

/* Variable Declarations */

/* Variable Definitions */
static SFc2_fcn2InstanceStruct chartInstance;

/* Function Declarations */
static void initialize_c2_fcn2(void);
static void initialize_params_c2_fcn2(void);
static void enable_c2_fcn2(void);
static void disable_c2_fcn2(void);
static void finalize_c2_fcn2(void);
static void sf_c2_fcn2(void);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber);
static void init_dsm_address_info(void);

/* Function Definitions */
static void initialize_c2_fcn2(void)
{
  _sfTime_ = (real_T)ssGetT(chartInstance.S);
}

static void initialize_params_c2_fcn2(void)
{
}

static void enable_c2_fcn2(void)
{
  _sfTime_ = (real_T)ssGetT(chartInstance.S);
}

static void disable_c2_fcn2(void)
{
  _sfTime_ = (real_T)ssGetT(chartInstance.S);
}

static void finalize_c2_fcn2(void)
{
}

static void sf_c2_fcn2(void)
{
  real_T *c2_y;
  real_T *c2_u;
  c2_y = (real_T *)ssGetOutputPortSignal(chartInstance.S, 1);
  c2_u = (real_T *)ssGetInputPortSignal(chartInstance.S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance.S);

  /*  This block supports the Embedded MATLAB subset. */
  /*  See the help menu for details.  */
  /*      y= u; */
  *c2_y = *c2_u;
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber)
{
}

const mxArray *sf_c2_fcn2_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo = NULL;
  c2_nameCaptureInfo = NULL;
  sf_mex_assign(&c2_nameCaptureInfo, sf_mex_create(NULL, "nameCaptureInfo", 0,
    0U, 1U, 2, 0, 1));
  return c2_nameCaptureInfo;
}

static void init_dsm_address_info(void)
{
}

/* SFunction Glue Code */
void sf_c2_fcn2_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2652502278U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(678159358U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1558451767U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(33149878U);
}

mxArray *sf_c2_fcn2_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,4,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateDoubleMatrix(4,1,mxREAL);
    double *pr = mxGetPr(mxChecksum);
    pr[0] = (double)(2094808043U);
    pr[1] = (double)(2147633222U);
    pr[2] = (double)(2948463122U);
    pr[3] = (double)(1580736031U);
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "aliasId", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,3,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"aliasId",mxCreateDoubleScalar(0));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "aliasId", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,3,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"aliasId",mxCreateDoubleScalar(0));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "aliasId", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,3,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"aliasId",mxCreateDoubleScalar(0));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  return(mxAutoinheritanceInfo);
}

static void sf_opaque_initialize_c2_fcn2(void *chartInstanceVar)
{
  initialize_params_c2_fcn2();
  initialize_c2_fcn2();
}

static void sf_opaque_enable_c2_fcn2(void *chartInstanceVar)
{
  enable_c2_fcn2();
}

static void sf_opaque_disable_c2_fcn2(void *chartInstanceVar)
{
  disable_c2_fcn2();
}

static void sf_opaque_gateway_c2_fcn2(void *chartInstanceVar)
{
  sf_c2_fcn2();
}

static void sf_opaque_terminate_c2_fcn2(void *chartInstanceVar)
{
  if (sim_mode_is_rtw_gen(chartInstance.S) || sim_mode_is_external
      (chartInstance.S)) {
    sf_clear_rtw_identifier(chartInstance.S);
  }

  finalize_c2_fcn2();
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_fcn2(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c2_fcn2();
  }
}

static void mdlSetWorkWidths_c2_fcn2(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable("fcn2","fcn2",2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop("fcn2","fcn2",2,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop("fcn2","fcn2",2,
      "gatewayCannotBeInlinedMultipleTimes"));
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,"fcn2","fcn2",2,2);
      sf_mark_chart_reusable_outputs(S,"fcn2","fcn2",2,1);
    }

    sf_set_rtw_dwork_info(S,"fcn2","fcn2",2);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
    ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  }

  ssSetChecksum0(S,(310977542U));
  ssSetChecksum1(S,(943118746U));
  ssSetChecksum2(S,(2557401175U));
  ssSetChecksum3(S,(731286671U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
}

static void mdlRTW_c2_fcn2(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    sf_write_symbol_mapping(S, "fcn2", "fcn2",2);
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_fcn2(SimStruct *S)
{
  chartInstance.chartInfo.chartInstance = NULL;
  chartInstance.chartInfo.isEMLChart = 1;
  chartInstance.chartInfo.chartInitialized = 0;
  chartInstance.chartInfo.sFunctionGateway = sf_opaque_gateway_c2_fcn2;
  chartInstance.chartInfo.initializeChart = sf_opaque_initialize_c2_fcn2;
  chartInstance.chartInfo.terminateChart = sf_opaque_terminate_c2_fcn2;
  chartInstance.chartInfo.enableChart = sf_opaque_enable_c2_fcn2;
  chartInstance.chartInfo.disableChart = sf_opaque_disable_c2_fcn2;
  chartInstance.chartInfo.zeroCrossings = NULL;
  chartInstance.chartInfo.outputs = NULL;
  chartInstance.chartInfo.derivatives = NULL;
  chartInstance.chartInfo.mdlRTW = mdlRTW_c2_fcn2;
  chartInstance.chartInfo.mdlStart = mdlStart_c2_fcn2;
  chartInstance.chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c2_fcn2;
  chartInstance.chartInfo.extModeExec = NULL;
  chartInstance.chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance.chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance.chartInfo.storeCurrentConfiguration = NULL;
  chartInstance.S = S;
  ssSetUserData(S,(void *)(&(chartInstance.chartInfo)));/* register the chart instance with simstruct */
  if (!sim_mode_is_rtw_gen(S)) {
    init_dsm_address_info();
  }
}

void c2_fcn2_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_fcn2(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_fcn2(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_fcn2(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_fcn2_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
