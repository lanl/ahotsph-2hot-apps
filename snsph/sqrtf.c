
float
sqrtf_fast(float x)
{
  float __value, __arg = (x);
  asm ("fsqrts %0,%1" : "=f" (__value) : "f" (__arg));
  return __value;
}


float
recip8bit(float x)
{
  float __value, __arg = (x);
  asm ("fres %0,%1" : "=f" (__value) : "f" (__arg));
  return __value;
}


/* Really 5bit */
float
recipsqrt8bit(float x)
{
  float __value, __arg = (x);
  asm ("frsqrte %0,%1" : "=f" (__value) : "f" (__arg));
  return __value;
}

