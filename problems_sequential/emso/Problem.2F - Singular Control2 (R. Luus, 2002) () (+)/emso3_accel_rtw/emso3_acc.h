/*
 * emso3_acc.h
 *
 * Real-Time Workshop code generation for Simulink model "emso3_acc.mdl".
 *
 * Model Version              : 1.124
 * Real-Time Workshop version : 7.1  (R2008a)  23-Jan-2008
 * C source code generated on : Fri Aug 16 16:24:56 2013
 */
#ifndef RTW_HEADER_emso3_acc_h_
#define RTW_HEADER_emso3_acc_h_
#ifndef emso3_acc_COMMON_INCLUDES_
# define emso3_acc_COMMON_INCLUDES_
#include <stdlib.h>
#define S_FUNCTION_NAME                simulink_only_sfcn
#define S_FUNCTION_LEVEL               2
#define RTW_GENERATED_S_FUNCTION
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_nonfinite.h"
#endif                                 /* emso3_acc_COMMON_INCLUDES_ */

#include "emso3_acc_types.h"

/* Block signals (auto storage) */
typedef struct {
  real_T B_0_0_0[3];                   /* '<Root>/EMSO model' */
} BlockIO;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  void *EMSOmodel_PWORK[13];           /* '<Root>/EMSO model' */
  struct {
    void *LoggedData;
  } Scope1_PWORK;                      /* '<Root>/Scope1' */

  struct {
    void *LoggedData;
  } Scope2_PWORK;                      /* '<Root>/Scope2' */

  struct {
    void *LoggedData;
  } Scope3_PWORK;                      /* '<Root>/Scope3' */

  struct {
    void *LoggedData;
  } Scope4_PWORK;                      /* '<Root>/Scope4' */

  int_T EMSOmodel_IWORK[4];            /* '<Root>/EMSO model' */
  char pad_EMSOmodel_IWORK[4];
} D_Work;

/* Continuous states (auto storage) */
typedef struct {
  real_T EMSOmodel_CSTATE[3];          /* '<Root>/EMSO model' */
} ContinuousStates;

/* State derivatives (auto storage) */
typedef struct {
  real_T EMSOmodel_CSTATE[3];          /* '<Root>/EMSO model' */
} StateDerivatives;

/* State disabled  */
typedef struct {
  boolean_T EMSOmodel_CSTATE[3];       /* '<Root>/EMSO model' */
} StateDisabled;

/* External inputs (root inport signals with auto storage) */
typedef struct {
  real_T In1[4];                       /* '<Root>/B_-1_-1' */
} ExternalUPtrs;

/* External outputs (root outports fed by signals with auto storage) */
typedef struct {
  real_T *B_0_1[3];                    /* '<Root>/B_0_1' */
} ExternalOutputs;

/* Parameters (auto storage) */
struct Parameters {
  real_T P_0[2];                       /* Computed Parameter: P1Size
                                        * '<Root>/EMSO model'
                                        */
  real_T P_1[14];                      /* Computed Parameter: P1
                                        * '<Root>/EMSO model'
                                        */
  real_T P_2[2];                       /* Computed Parameter: P2Size
                                        * '<Root>/EMSO model'
                                        */
  real_T P_3;                          /* Expression: 4
                                        * '<Root>/EMSO model'
                                        */
  real_T P_4[2];                       /* Computed Parameter: P3Size
                                        * '<Root>/EMSO model'
                                        */
  real_T P_5;                          /* Expression: 3
                                        * '<Root>/EMSO model'
                                        */
};

extern Parameters rtDefaultParameters; /* parameters */

#endif                                 /* RTW_HEADER_emso3_acc_h_ */
