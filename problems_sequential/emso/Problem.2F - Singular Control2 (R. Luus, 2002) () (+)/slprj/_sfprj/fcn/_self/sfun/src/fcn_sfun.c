/* Include files */

#include "fcn_sfun.h"
#include "c2_fcn.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint8_T _sfEvent_;
uint32_T _fcnMachineNumber_;
real_T _sfTime_;

/* Function Declarations */

/* Function Definitions */
void fcn_initializer(void)
{
  _sfEvent_ = CALL_EVENT;
}

void fcn_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_fcn_method_dispatcher(SimStruct *simstructPtr, unsigned int
  chartFileNumber, int_T method, void *data)
{
  if (chartFileNumber==2) {
    c2_fcn_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_fcn_process_check_sum_call( int nlhs, mxArray * plhs[], int nrhs,
  const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1205617842U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3808315393U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3958808787U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(816851450U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3660826754U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3611417902U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2806998225U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(875321426U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 2:
        {
          extern void sf_c2_fcn_get_check_sum(mxArray *plhs[]);
          sf_c2_fcn_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(64876186U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1217120569U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(4094067315U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1782201511U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1868985424U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(800849192U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(339174905U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1968041558U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_fcn_autoinheritance_info( int nlhs, mxArray * plhs[], int nrhs,
  const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  if (nrhs<2 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 2:
      {
        extern mxArray *sf_c2_fcn_get_autoinheritance_info(void);
        plhs[0] = sf_c2_fcn_get_autoinheritance_info();
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_fcn_get_eml_resolved_functions_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 2:
      {
        extern const mxArray *sf_c2_fcn_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_fcn_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

void fcn_debug_initialize(void)
{
  _fcnMachineNumber_ = sf_debug_initialize_machine("fcn","sfun",0,1,0,0,0);
  sf_debug_set_machine_event_thresholds(_fcnMachineNumber_,0,0);
  sf_debug_set_machine_data_thresholds(_fcnMachineNumber_,0);
}

void fcn_register_exported_symbols(SimStruct* S)
{
}
