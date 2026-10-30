/*
 * emso3_acc_private.h
 *
 * Real-Time Workshop code generation for Simulink model "emso3_acc.mdl".
 *
 * Model Version              : 1.124
 * Real-Time Workshop version : 7.1  (R2008a)  23-Jan-2008
 * C source code generated on : Fri Aug 16 16:24:56 2013
 */
#ifndef RTW_HEADER_emso3_acc_private_h_
#define RTW_HEADER_emso3_acc_private_h_
#include "rtwtypes.h"
#ifndef RTW_COMMON_DEFINES_
# define RTW_COMMON_DEFINES_
#define rt_VALIDATE_MEMORY(S, ptr)     if(!(ptr)) {\
 ssSetErrorStatus(S, RT_MEMORY_ALLOCATION_ERROR);\
 return;\
 }
#if !defined(_WIN32)
#define rt_FREE(ptr)                   if((ptr) != NULL) {\
 free((ptr));\
 (ptr) = NULL;\
 }
#else

/* Visual and other windows compilers delcare free without const! */
#define rt_FREE(ptr)                   if((ptr) != NULL) {\
 free((void *)(ptr));\
 (ptr) = NULL;\
 }
#endif
#endif                                 /* RTW_COMMON_DEFINES_ */

#ifndef __RTWTYPES_H__
#error This file requires rtwtypes.h to be included
#else
#ifdef TMWTYPES_PREVIOUSLY_INCLUDED
#error This file requires rtwtypes.h to be included before tmwtypes.h
#else

/* Check for inclusion of an incorrect version of rtwtypes.h */
#ifndef RTWTYPES_ID_C08S16I32L32N32F1
#error This code was generated with a different "rtwtypes.h" than the file included
#endif                                 /* RTWTYPES_ID_C08S16I32L32N32F1 */
#endif                                 /* TMWTYPES_PREVIOUSLY_INCLUDED */
#endif                                 /* __RTWTYPES_H__ */
#endif                                 /* RTW_HEADER_emso3_acc_private_h_ */
