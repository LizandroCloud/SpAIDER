#include "rt_nonfinite.h"
#define NumBitsPerChar                 8

real_T rtInf;
real_T rtMinusInf;
real_T rtNaN;
real32_T rtInfF;
real32_T rtMinusInfF;
real32_T rtNaNF;
void rt_InitInfAndNaN(size_t realSize)
{
  size_t bitsPerReal = realSize * (NumBitsPerChar);
  uint16_T one = 1;
  enum {
    LittleEndian,
    BigEndian
  } machByteOrder = (*((uint8_T *) &one) == 1) ? LittleEndian : BigEndian;
  switch (machByteOrder) {
   case LittleEndian:
    {
      typedef struct {
        uint32_T fraction : 23;
        uint32_T exponent : 8;
        uint32_T sign : 1;
      } LittleEndianIEEESingle;

      union {
        LittleEndianIEEESingle bitVal;
        real32_T fltVal;
      } tmpValF;

      typedef struct {
        struct {
          uint32_T fraction2;
        } wordH;

        struct {
          uint32_T fraction1 : 20;
          uint32_T exponent : 11;
          uint32_T sign : 1;
        } wordL;
      } LittleEndianIEEEDouble;

      union {
        LittleEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpValF.bitVal.sign = 0;
      tmpValF.bitVal.exponent = 0xFF;
      tmpValF.bitVal.fraction = 0;
      rtInfF = tmpValF.fltVal;
      tmpValF.bitVal.sign = 1;
      rtMinusInfF = tmpValF.fltVal;
      tmpValF.bitVal.fraction = 0x400000;
      rtNaNF = tmpValF.fltVal;
      if (bitsPerReal == 32) {
        rtInf = rtInfF;
        rtMinusInf = rtMinusInfF;
        rtNaN = rtNaNF;
      } else {
        tmpVal.bitVal.wordL.sign = 0;
        tmpVal.bitVal.wordL.exponent = 0x7FF;
        tmpVal.bitVal.wordL.fraction1 = 0;
        tmpVal.bitVal.wordH.fraction2 = 0;
        rtInf = tmpVal.fltVal;
        tmpVal.bitVal.wordL.sign = 1;
        rtMinusInf = tmpVal.fltVal;
        tmpVal.bitVal.wordL.fraction1 = 0x80000;
        rtNaN = tmpVal.fltVal;
      }
    }
    break;

   case BigEndian:
    {
      typedef struct {
        uint32_T sign : 1;
        uint32_T exponent : 8;
        uint32_T fraction : 23;
      } BigEndianIEEESingle;

      union {
        BigEndianIEEESingle bitVal;
        real32_T fltVal;
      } tmpValF;

      typedef struct {
        struct {
          uint32_T sign : 1;
          uint32_T exponent : 11;
          uint32_T fraction1 : 20;
        } wordL;

        struct {
          uint32_T fraction2;
        } wordH;
      } BigEndianIEEEDouble;

      union {
        BigEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpValF.bitVal.sign = 0;
      tmpValF.bitVal.exponent = 0xFF;
      tmpValF.bitVal.fraction = 0;
      rtInfF = tmpValF.fltVal;
      tmpValF.bitVal.sign = 1;
      rtMinusInfF = tmpValF.fltVal;
      tmpValF.bitVal.fraction = 0x400000;
      rtNaNF = tmpValF.fltVal;
      tmpValF.bitVal.sign = 0;
      tmpValF.bitVal.fraction = 0x7fffff;
      rtNaNF = tmpValF.fltVal;
      if (bitsPerReal == 32) {
        rtInf = rtInfF;
        rtMinusInf = rtMinusInfF;
        rtNaN = rtNaNF;
      } else {
        tmpVal.bitVal.wordL.sign = 0;
        tmpVal.bitVal.wordL.exponent = 0x7FF;
        tmpVal.bitVal.wordL.fraction1 = 0;
        tmpVal.bitVal.wordH.fraction2 = 0;
        rtInf = tmpVal.fltVal;
        tmpVal.bitVal.wordL.sign = 1;
        rtMinusInf = tmpVal.fltVal;
        tmpVal.bitVal.wordL.sign = 0;
        tmpVal.bitVal.wordL.fraction1 = 0xFFFFF;
        tmpVal.bitVal.wordH.fraction2 = 0xFFFFFFFF;
        rtNaN = tmpVal.fltVal;
      }
    }
    break;
  }
}

boolean_T rtIsInf(real_T value)
{
  return ((value==rtInf || value==rtMinusInf) ? 1U : 0U);
}

boolean_T rtIsInfF(real32_T value)
{
  return(((value)==rtInfF || (value)==rtMinusInfF) ? 1U : 0U);
}

boolean_T rtIsNaN(real_T value)
{
  return _isnan(value)? 1U:0U;
}

boolean_T rtIsNaNF(real32_T value)
{
  return _isnan((real_T)value)? 1U:0U;
}

#undef NumBitsPerChar

