# 1 "do_grav_fast.c"
# 32 "do_grav_fast.c"
// Things to do:
//  Check doubleword alignments
//  Avoid need for padding
# 36 "do_grav_fast.c"
// parameters
// const pointers r18 and r22 are free to use after they have been dereferenced
# 50 "do_grav_fast.c"
// The list that is processed needs padding at the r17 to account for
// r17 - pp not divisible by twelve, and also to provide valid data that
// is read on the last iteration, but not processed.
// This padding may be provided by 20 zeros at the r17 of the list.
//
// void do_grav(float *r16, float *r17, const float *r18, float *r19,
// 	     float *r20, float *r21, const float *r22, int *r23)
# 58 "do_grav_fast.c"
	.align	32
_do_grav_fast::
.a4 = 0
.f4 = 32
	addu -(.a4+.f4), sp, sp
	st.l r1, 0(r2)
# 65 "do_grav_fast.c"
	subu r16, r17, r0		// pp - r17
	bc .do_nothing
# 68 "do_grav_fast.c"
	fst.d f2, 8(r2)		// store f2-f7, left unchanged by subroutine
	fst.d f4, 16(r2)
	fst.d f6, 24(r2)
# 72 "do_grav_fast.c"
	// r16, float *r16
	// r17, float *r17
	// r18, const float *r18
	// r19, float *r19
	// r20, float *r20
	// r21, float *r21
	// r22, const float *r22
	// r23, int *r23
	// r24, int *tot_interact
# 82 "do_grav_fast.c"
	fld.l r0(r22), f8
					orh 0x3f00,r0,r31	// 0.5(float)
	fld.l 0(r18), f2
					ixfr r31,f31
	fld.l 4(r18), f3
					orh 0x3fc0,r0,r31	// 1.5(float)
	fld.l 8(r18), f4
					ixfr r31,f30
	fld.l 0(r20), f5
					addu 16, r17, r17
	fld.l 4(r20), f6
	fld.l 8(r20), f7
	fld.l r0(r21), f9
	fld.l r0(r19), f22
	ld.l r0(r23), r18			// use r18
# 98 "do_grav_fast.c"
	pfld.d 0(r16), f0		 	// no autoincrement for first
	pfld.d 8(r16)++, f0
	pfld.d 8(r16)++, f0
# 102 "do_grav_fast.c"
	pfld.d 8(r16)++, f10
	pfld.d 8(r16)++, f12
	pfld.d 8(r16)++, f14
	pfld.d 8(r16)++, f16
	pfld.d 8(r16)++, f18
	pfld.d 8(r16)++, f20
# 109 "do_grav_fast.c"
// 106 clocks in main loop
.main_loop:
	pfsub.ss f11, f2, f0
	pfsub.ss f15, f2, f0
	pfsub.ss f19, f2, f0
	pfsub.ss f12, f3, f11
	pfsub.ss f16, f3, f15
	pfsub.ss f20, f3, f19
	pfsub.ss f13, f4, f12
	pfsub.ss f17, f4, f16
	pfsub.ss f21, f4, f20
	pfmul.ss f11, f11, f0
	pfmul.ss f15, f15, f0
	pfmul.ss f19, f19, f0
	pfmul.ss f12, f12, f27
	pfmul.ss f16, f16, f28
	pfmul.ss f20, f20, f29
	pfadd.ss f27, f8, f13
	pfadd.ss f28, f8, f17
	pfadd.ss f29, f8, f21
	m12apm.ss f13, f13, f0			// adder + multiplier
	m12apm.ss f17, f17, f0
	m12apm.ss f21, f21, f0
	m12apm.ss f0, f0, f0			// adder + multiplier
	m12apm.ss f0, f0, f0
	m12apm.ss f0, f0, f0
# 136 "do_grav_fast.c"
	r2apt.ss f0, f0, f27			// clear KR for later
	r2apt.ss f0, f0, f28
	r2apt.ss f0, f0, f29
# 140 "do_grav_fast.c"
				// frsqr screws up the M pipe in sim860
				// pipelined sqrtf3, 39 clocks
	frsqr.ss f27, f23		// initial guess of inverse sqrt
	frsqr.ss f28, f24
	frsqr.ss f29, f25
	pfmul.ss f31, f27, f0	// 0.5 * x
	pfmul.ss f31, f28, f0
	pfmul.ss f31, f29, f0
	pfmul.ss f23, f23, f27	// yold * yold
	pfmul.ss f24, f24, f28
	pfmul.ss f25, f25, f29
	pfmul.ss f27, f26, f26	// (0.5 * x) * (yold * yold)
	pfmul.ss f28, f26, f26
	pfmul.ss f29, f26, f26
	mr2s1.ss  f30, f0, f0	// 1.5 - ( ) // subtract src1, M result
	mr2s1.ss  f30, f0, f0
	mr2s1.ss  f30, f0, f0
	pfadd.ss f0, f0, f26	// no other way to get adder into M
	pfmul.ss f23, f26, f0	// yold * ( )
	pfadd.ss f0, f0, f26
	pfmul.ss f24, f26, f0
	pfadd.ss f0, f0, f26
	pfmul.ss f25, f26, f0
	pfmul.ss f27, f23, f23	// (0.5 * x) * yold
	pfmul.ss f28, f24, f24
	pfmul.ss f29, f25, f25
	pfmul.ss f23, f26, f26	// yold * ( )
	pfmul.ss f24, f26, f26
	pfmul.ss f25, f26, f26
	mr2s1.ss  f30, f0, f0	// 1.5 - ( ) // subtract src1, M result
	mr2s1.ss  f30, f0, f0
	mr2s1.ss  f30, f0, f0
	pfadd.ss f0, f0, f26
	pfmul.ss f23, f26, f0	// yold * ( )
	pfadd.ss f0, f0, f26
	pfmul.ss f24, f26, f0
	pfadd.ss f0, f0, f26
	pfmul.ss f25, f26, f0
	pfmul.ss f0, f0, f27	// result
	pfmul.ss f0, f0, f28
	pfmul.ss f0, f0, f29
# 182 "do_grav_fast.c"
		// pflds should not be executed on the last cycle
		// need to pad list, or write more code to handle the r17
	pfmul.ss f27, f27, f0
	pfmul.ss f28, f28, f0
	pfmul.ss f29, f29, f0
	pfmul.ss f27, f10, f23
	pfmul.ss f28, f14, f24
	pfmul.ss f29, f18, f25
	pfmul.ss f23, f27, f27
	pfmul.ss f24, f28, f28
	pfmul.ss f25, f29, f29
	pfadd.ss f10, f22, f0
	pfsub.ss f9, f27,  f0
	pfadd.ss f0, f0, f0
	pfadd.ss f14, f22, f22
	pfsub.ss f9, f28,  f9
	pfadd.ss f0, f0, f0
	pfadd.ss f18, f22, f22
	pfsub.ss f9, f29,  f9
	pfmul.ss f11, f23, f23
	pfmul.ss f12, f23, f24
	pfmul.ss f13, f23, f25
	r2p1.ss  f5, f0, f0		// src1 + M result, advance adder
	r2p1.ss  f6, f0, f22
	r2p1.ss  f7, f0, f9
.dual
	pfmul.ss f15, f24, f0
	pfmul.ss f16, f24, f0
	pfmul.ss f17, f24, f0
					pfld.d 8(r16)++, f10
	m12apm.ss  f0, f0, f0		// A result + M result, advance adder
					subu r16, r17, r0		// r16 - r17
	m12apm.ss  f0, f0, f0
					pfld.d 8(r16)++, f12
	m12apm.ss  f0, f0, f0	
					nop
	pfmul.ss f19, f25, f0
					pfld.d 8(r16)++, f14
	pfmul.ss f20, f25, f0
					nop
	pfmul.ss f21, f25, f0
					pfld.d 8(r16)++, f16
	m12apm.ss  f0, f0, f0		// A result + M result, advance adder
					nop
	m12apm.ss  f0, f0, f0
					pfld.d 8(r16)++, f18
.enddual
	m12apm.ss  f0, f0, f0	
					nop
	pfadd.ss f0, f0, f5
					pfld.d 8(r16)++, f20
	pfadd.ss f0, f0, f6
	pfadd.ss f0, f0, f7
# 236 "do_grav_fast.c"
	bnc.t .main_loop
	nop
# 239 "do_grav_fast.c"
	fst.l f5, 0(r20)	// f5
	fst.l f6, 4(r20)	// f6
	fst.l f7, 8(r20)	// f7
	fst.l f9, r0(r21)	// f9
	fst.l f22, r0(r19)	// f22
	st.l r18, 0(r23)
# 246 "do_grav_fast.c"
	fld.d 8(r2), f2
	fld.d 16(r2), f4
	fld.d 24(r2), f6
# 250 "do_grav_fast.c"
.do_nothing:
	ld.l 0(r2), r1
	bri r1
	addu (.a4+.f4), sp, sp
