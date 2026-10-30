/*
 * This file is not available for use in any application other than as a
 * MATLAB(R) MEX file for use with the Simulink(R) product.
 */

/*
 * emso3_acc.c
 *
 * Real-Time Workshop code generation for Simulink model "emso3_acc.mdl".
 *
 * Model Version              : 1.124
 * Real-Time Workshop version : 7.1  (R2008a)  23-Jan-2008
 * C source code generated on : Fri Aug 16 16:24:56 2013
 */
#include <math.h>
#include "emso3_acc.h"
#include "emso3_acc_private.h"
#include <stdio.h>
#include "simstruc.h"
#include "fixedpoint.h"
#define CodeFormat                     S-Function
#define AccDefine1                     Accelerator_S-Function

/* Outputs for root system: '<Root>' */
static void mdlOutputs(SimStruct *S, int_T tid)
{
  /* Level2 S-Function Block: '<Root>/B_0_0' (emso_sf_c) */
  /* Call into Simulink for MEX-version of S-function */
  ssCallAccelRunBlock(S, 0, 0, SS_CALL_MDL_OUTPUTS);
  ssCallAccelRunBlock(S, 0, 1, SS_CALL_MDL_OUTPUTS);

  /* Scope: '<Root>/Scope1' */

  /* Call into Simulink for Scope */
  ssCallAccelRunBlock(S, 0, 2, SS_CALL_MDL_OUTPUTS);

  /* Scope: '<Root>/Scope2' */

  /* Call into Simulink for Scope */
  ssCallAccelRunBlock(S, 0, 3, SS_CALL_MDL_OUTPUTS);

  /* Scope: '<Root>/Scope3' */

  /* Call into Simulink for Scope */
  ssCallAccelRunBlock(S, 0, 4, SS_CALL_MDL_OUTPUTS);

  /* Scope: '<Root>/Scope4' */

  /* Call into Simulink for Scope */
  ssCallAccelRunBlock(S, 0, 5, SS_CALL_MDL_OUTPUTS);

  /* tid is required for a uniform function interface. This system
   * is single rate, and in this case, tid is not accessed. */
  UNUSED_PARAMETER(tid);
}

/* Update for root system: '<Root>' */
#define MDL_UPDATE

static void mdlUpdate(SimStruct *S, int_T tid)
{
  /* tid is required for a uniform function interface. This system
   * is single rate, and in this case, tid is not accessed. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
#define MDL_DERIVATIVES

static void mdlDerivatives(SimStruct *S)
{
  /* Level2 S-Function Block: '<Root>/B_0_0' (emso_sf_c) */
  /* Call into Simulink for MEX-version of S-function */
  ssCallAccelRunBlock(S, 0, 0, SS_CALL_MDL_DERIVATIVES);
}

/* Function to initialize sizes */
static void mdlInitializeSizes(SimStruct *S)
{
  /* checksum */
  ssSetChecksumVal(S, 0, 2335424531U);
  ssSetChecksumVal(S, 1, 2562509485U);
  ssSetChecksumVal(S, 2, 3874560507U);
  ssSetChecksumVal(S, 3, 4218365658U);

  /* options */
  ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);

  /* Accelerator check memory map size match for DWork */
  if (ssGetSizeofDWork(S) != sizeof(D_Work)) {
    ssSetErrorStatus(S,"Unexpected error: Internal DWork sizes do "
                     "not match for accelerator mex file.");
  }

  /* Accelerator check memory map size match for BlockIO */
  if (ssGetSizeofGlobalBlockIO(S) != sizeof(BlockIO)) {
    ssSetErrorStatus(S,"Unexpected error: Internal BlockIO sizes do "
                     "not match for accelerator mex file.");
  }

  /* model parameters */
  _ssSetDefaultParam(S, (real_T *) &rtDefaultParameters);

  /* non-finites */
  rt_InitInfAndNaN(sizeof(real_T));
}

/* Empty mdlInitializeSampleTimes function (never called) */
static void mdlInitializeSampleTimes(SimStruct *S)
{
}

/* Empty mdlTerminate function (never called) */
static void mdlTerminate(SimStruct *S)
{
}

/* MATLAB MEX Glue */
#include "simulink.c"
