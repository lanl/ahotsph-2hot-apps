// This is the C callable assembler code for sqrtf function 
// This program is written by Edith Huang at Caltech/JPL, Sept.,1991 
// This one takes 1.21musec (johns)
// uses registers f8, f17, f18, f19, f20, f21, f22, f24
// f8    // Contains X
_sqrtf_fast::

// Test to see if X>0
     
     pfgt.ss f8, f0, f0    // if X > 0

     bc.t   .BINSQRT
     orh    0x3f00,r0,r31      // 0.5(float) -> r31

     pfeq.ss f8,f0,f0     // if x=0
     bc   .BZERO

     orh    0x7fc0,r0,r31
     ixfr   r31,f8          // return signalling NaN for error

     bri r1
     nop

.BZERO:

     fmov.ss f0,f8      //  return 0 when x=0
     bri r1
     nop

.BINSQRT:
     frsqr.ss  f8, f18     // initial guess of inverse sqrt
     ixfr   r31,f17		// 0.5 -> f17
     fmul.ss   f18, f18, f19     // yold **2 -> f19
     orh    0x3fc0,r0,r31         // 1.5(float) -> r31
// changing f17,f8 to f8,f17 on the next line slows things down by 5%!?
     fmul.ss   f17,f8, f20      // 0.5*x -> f20
     ixfr   r31, f22		// 1.5 -> f22
     fmul.ss   f19,f20 , f21     // 0.5*x*yold*yold -> f21
     fsub.ss   f22, f21, f24      // 1.5 - 0.5 *x*yold**2 
     fmul.ss   f18,f24,f18       //  yold(1.5-0.5*x*yold**2)

// 2nd iteration

     fmul.ss   f18,f18,f19	// yold**2
     fmul.ss   f20,f19 , f21     // 0.5*x*yold*yold
     fsub.ss   f22, f21, f24      // 1.5 - 0.5 *x*yold**2 
     fmul.ss   f18,f24,f18       //  yold(1.5-0.5*x*yold**2)
     fmul.ss   f18,f8,f8

// result

     bri  r1
     nop

