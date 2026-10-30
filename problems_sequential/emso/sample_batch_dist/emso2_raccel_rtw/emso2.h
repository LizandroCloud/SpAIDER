#ifndef RTW_HEADER_emso2_h_
#define RTW_HEADER_emso2_h_
#ifndef emso2_COMMON_INCLUDES_
# define emso2_COMMON_INCLUDES_
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "raccel.h"
#include "rt_logging.h"
#include "dt_info.h"
#include "ext_work.h"
#include "rt_nonfinite.h"
#endif

#include "emso2_types.h"
#define MODEL_NAME                     emso2
#define NSAMPLE_TIMES                  (1)
#define NINPUTS                        (1)
#define NOUTPUTS                       (4)
#define NBLOCKIO                       (1)
#define NUM_ZC_EVENTS                  (0)
#ifndef NCSTATES
# define NCSTATES                      (0)
#elif NCSTATES != 0
# error Invalid specification of NCSTATES defined in compiler command
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        NULL
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)
#endif

typedef struct {
  real_T EMSOmodel[4];
} BlockIO;

typedef struct {
  void *EMSOmodel_PWORK[13];
  struct {
    void *LoggedData;
  } Scope_PWORK;

  struct {
    void *LoggedData;
  } Scope2_PWORK;

  struct {
    void *LoggedData;
  } Scope3_PWORK;

  struct {
    void *LoggedData;
  } Scope4_PWORK;

  struct {
    void *LoggedData;
  } Scope1_PWORK;

  int_T EMSOmodel_IWORK[3];
} D_Work;

typedef struct {
  real_T In1;
} ExternalInputs;

typedef struct {
  real_T Out1[4];
} ExternalOutputs;

struct Parameters {
  real_T Constant3_Value;
  real_T EMSOmodel_P1_Size[2];
  real_T EMSOmodel_P1[14];
  real_T EMSOmodel_P2_Size[2];
  real_T EMSOmodel_P2;
  real_T EMSOmodel_P3_Size[2];
  real_T EMSOmodel_P3;
  real_T EMSOmodel_P4_Size[2];
  real_T EMSOmodel_P4;
};

extern Parameters rtP;
extern const char *RT_MEMORY_ALLOCATION_ERROR;
extern BlockIO rtB;
extern D_Work rtDWork;
extern ExternalInputs rtU;
extern ExternalOutputs rtY;
extern SimStruct *const rtS;
extern const int_T gblNumToFiles;
extern const int_T gblNumFrFiles;
extern const int_T gblNumFrWksBlocks;
extern rtInportTUtable *gblInportTUtables;
extern const char *gblInportFileName;
extern const int_T gblNumRootInportBlks;
extern const int_T gblNumModelInputs;
extern const int_T gblInportDataTypeIdx[];
extern const int_T gblInportDims[];
extern const int_T gblInportComplex[];
extern const int_T gblInportInterpoFlag[];
extern const int_T gblInportContinuous[];

#endif
