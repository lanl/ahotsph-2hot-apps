#define ppos0 f2
#define ppos1 f3
#define ppos2 f4
#define a0 f5
#define a1 f6
#define a2 f7
#define eps2 f8
#define phi f9
#define m0 f10
#define x0 f11
#define y0 f12
#define z0 f13
#define m1 f14
#define x1 f15
#define y1 f16
#define z1 f17
#define m2 f18
#define x2 f19
#define y2 f20
#define z2 f21
#define mass f22
#define t0 f23
#define t1 f24
#define t2 f25
#define t3 f26
#define g0 f27
#define g1 f28
#define g2 f29
#define c1 f30
#define c2 f31

// Things to do:
//  Check doubleword alignments
//  Avoid need for padding

// parameters
// const pointers r18 and r22 are free to use after they have been dereferenced

#define p r16
#define end r17
#define pos0 r18
#define mass0 r19
#define acc0 r20
#define phi0 r21
#define eps2p r22
#define ncut r23

#define cut r18

// The list that is processed needs padding at the end to account for
// end - pp not divisible by twelve, and also to provide valid data that
// is read on the last iteration, but not processed.
// This padding may be provided by 20 zeros at the end of the list.
//
// void do_grav(float *p, float *end, const float *pos0, float *mass0,
// 	     float *acc0, float *phi0, const float *eps2p, int *ncut)

	.align	32
_do_grav_fast::
.a4 = 0
.f4 = 32
	addu -(.a4+.f4), sp, sp
	st.l r1, 0(r2)

	subu p, end, r0		// pp - end
	bc .do_nothing

	fst.d f2, 8(r2)		// store f2-f7, left unchanged by subroutine
	fst.d f4, 16(r2)
	fst.d f6, 24(r2)

	// r16, float *p
	// r17, float *end
	// r18, const float *pos0
	// r19, float *mass0
	// r20, float *acc0
	// r21, float *phi0
	// r22, const float *eps2p
	// r23, int *ncut
	// r24, int *tot_interact

	fld.l r0(eps2p), eps2
					orh 0x3f00,r0,r31	// 0.5(float)
	fld.l 0(pos0), ppos0
					ixfr r31,c2
	fld.l 4(pos0), ppos1
					orh 0x3fc0,r0,r31	// 1.5(float)
	fld.l 8(pos0), ppos2
					ixfr r31,c1
	fld.l 0(acc0), a0
					addu 16, end, end
	fld.l 4(acc0), a1
	fld.l 8(acc0), a2
	fld.l r0(phi0), phi
	fld.l r0(mass0), mass
	ld.l r0(ncut), cut			// use r18

	pfld.d 0(p), f0		 	// no autoincrement for first
	pfld.d 8(p)++, f0
	pfld.d 8(p)++, f0

	pfld.d 8(p)++, m0
	pfld.d 8(p)++, y0
	pfld.d 8(p)++, m1
	pfld.d 8(p)++, y1
	pfld.d 8(p)++, m2
	pfld.d 8(p)++, y2

// 106 clocks in main loop
.main_loop:
	pfsub.ss x0, ppos0, f0
	pfsub.ss x1, ppos0, f0
	pfsub.ss x2, ppos0, f0
	pfsub.ss y0, ppos1, x0
	pfsub.ss y1, ppos1, x1
	pfsub.ss y2, ppos1, x2
	pfsub.ss z0, ppos2, y0
	pfsub.ss z1, ppos2, y1
	pfsub.ss z2, ppos2, y2
	pfmul.ss x0, x0, f0
	pfmul.ss x1, x1, f0
	pfmul.ss x2, x2, f0
	pfmul.ss y0, y0, g0
	pfmul.ss y1, y1, g1
	pfmul.ss y2, y2, g2
	pfadd.ss g0, eps2, z0
	pfadd.ss g1, eps2, z1
	pfadd.ss g2, eps2, z2
	m12apm.ss z0, z0, f0			// adder + multiplier
	m12apm.ss z1, z1, f0
	m12apm.ss z2, z2, f0
	m12apm.ss f0, f0, f0			// adder + multiplier
	m12apm.ss f0, f0, f0
	m12apm.ss f0, f0, f0

	r2apt.ss f0, f0, g0			// clear KR for later
	r2apt.ss f0, f0, g1
	r2apt.ss f0, f0, g2

				// frsqr screws up the M pipe in sim860
				// pipelined sqrtf3, 39 clocks
	frsqr.ss g0, t0		// initial guess of inverse sqrt
	frsqr.ss g1, t1
	frsqr.ss g2, t2
	pfmul.ss c2, g0, f0	// 0.5 * x
	pfmul.ss c2, g1, f0
	pfmul.ss c2, g2, f0
	pfmul.ss t0, t0, g0	// yold * yold
	pfmul.ss t1, t1, g1
	pfmul.ss t2, t2, g2
	pfmul.ss g0, t3, t3	// (0.5 * x) * (yold * yold)
	pfmul.ss g1, t3, t3
	pfmul.ss g2, t3, t3
	mr2s1.ss  c1, f0, f0	// 1.5 - ( ) // subtract src1, M result
	mr2s1.ss  c1, f0, f0
	mr2s1.ss  c1, f0, f0
	pfadd.ss f0, f0, t3	// no other way to get adder into M
	pfmul.ss t0, t3, f0	// yold * ( )
	pfadd.ss f0, f0, t3
	pfmul.ss t1, t3, f0
	pfadd.ss f0, f0, t3
	pfmul.ss t2, t3, f0
	pfmul.ss g0, t0, t0	// (0.5 * x) * yold
	pfmul.ss g1, t1, t1
	pfmul.ss g2, t2, t2
	pfmul.ss t0, t3, t3	// yold * ( )
	pfmul.ss t1, t3, t3
	pfmul.ss t2, t3, t3
	mr2s1.ss  c1, f0, f0	// 1.5 - ( ) // subtract src1, M result
	mr2s1.ss  c1, f0, f0
	mr2s1.ss  c1, f0, f0
	pfadd.ss f0, f0, t3
	pfmul.ss t0, t3, f0	// yold * ( )
	pfadd.ss f0, f0, t3
	pfmul.ss t1, t3, f0
	pfadd.ss f0, f0, t3
	pfmul.ss t2, t3, f0
	pfmul.ss f0, f0, g0	// result
	pfmul.ss f0, f0, g1
	pfmul.ss f0, f0, g2

		// pflds should not be executed on the last cycle
		// need to pad list, or write more code to handle the end
	pfmul.ss g0, g0, f0
	pfmul.ss g1, g1, f0
	pfmul.ss g2, g2, f0
	pfmul.ss g0, m0, t0
	pfmul.ss g1, m1, t1
	pfmul.ss g2, m2, t2
	pfmul.ss t0, g0, g0
	pfmul.ss t1, g1, g1
	pfmul.ss t2, g2, g2
	pfadd.ss m0, mass, f0
	pfsub.ss phi, g0,  f0
	pfadd.ss f0, f0, f0
	pfadd.ss m1, mass, mass
	pfsub.ss phi, g1,  phi
	pfadd.ss f0, f0, f0
	pfadd.ss m2, mass, mass
	pfsub.ss phi, g2,  phi
	pfmul.ss x0, t0, t0
	pfmul.ss y0, t0, t1
	pfmul.ss z0, t0, t2
	r2p1.ss  a0, f0, f0		// src1 + M result, advance adder
	r2p1.ss  a1, f0, mass
	r2p1.ss  a2, f0, phi
.dual
	pfmul.ss x1, t1, f0
	pfmul.ss y1, t1, f0
	pfmul.ss z1, t1, f0
					pfld.d 8(p)++, m0
	m12apm.ss  f0, f0, f0		// A result + M result, advance adder
					subu p, end, r0		// p - end
	m12apm.ss  f0, f0, f0
					pfld.d 8(p)++, y0
	m12apm.ss  f0, f0, f0	
					nop
	pfmul.ss x2, t2, f0
					pfld.d 8(p)++, m1
	pfmul.ss y2, t2, f0
					nop
	pfmul.ss z2, t2, f0
					pfld.d 8(p)++, y1
	m12apm.ss  f0, f0, f0		// A result + M result, advance adder
					nop
	m12apm.ss  f0, f0, f0
					pfld.d 8(p)++, m2
.enddual
	m12apm.ss  f0, f0, f0	
					nop
	pfadd.ss f0, f0, a0
					pfld.d 8(p)++, y2
	pfadd.ss f0, f0, a1
	pfadd.ss f0, f0, a2

	bnc.t .main_loop
	nop

	fst.l a0, 0(acc0)	// a0
	fst.l a1, 4(acc0)	// a1
	fst.l a2, 8(acc0)	// a2
	fst.l phi, r0(phi0)	// phi
	fst.l mass, r0(mass0)	// mass
	st.l cut, 0(ncut)

	fld.d 8(r2), f2
	fld.d 16(r2), f4
	fld.d 24(r2), f6

.do_nothing:
	ld.l 0(r2), r1
	bri r1
	addu (.a4+.f4), sp, sp
