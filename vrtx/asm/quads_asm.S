/* Xd and tb share regs */
#define Xd0 f2
#define Xd1 f3
#define Xd2 f4
#define tb0 f2
#define tb1 f3
#define tb2 f4

#define cf0 f5
#define XQ0 f5
#define XQ1 f6
#define XQ2 f7
#define wz0 f5
#define wz1 f6
#define wz2 f7

#define x f8
#define y f9
#define z f10

#define Qz0 f11
#define Qz1 f15
#define Qz2 f20
#define xy f11
#define xz f15
#define yz f20
#define xzp0 f11
#define xzp1 f15
#define xzp2 f20

#define Qx0 f16
#define Qx1 f17
#define Qx2 f18
#define xx f16
#define yy f17
#define zz f18
#define xxp1 f16
#define yyp2 f17
#define zzp0 f18

#define Qy0 f19
#define Qy1 f22
#define Qy2 f23
#define m2x f19
#define m2y f22
#define m2z f23
#define xyp0 f19
#define xyp1 f22
#define xyp2 f23

#define cf1 f29
#define cf3 f27
#define cf5 f28
#define yzp0 f29
#define yzp1 f27
#define yzp2 f28

#define ta0 f12
#define ta1 f13
#define ta2 f14
#define st0 f12
#define st1 f13
#define st2 f14
#define xxp2 f12
#define yyp0 f13
#define zzp1 f14

#define t0 f30
#define t1 f21
#define t2 f31

#define cf3_2 f24
#define cf5_2 f25
#define cf7_2 f26
#define s0 f24
#define s1 f25
#define s2 f26
#define wy0 f24
#define wy1 f25
#define wy2 f26

#define St0 12(r16)
#define St1 16(r16)
#define St2 20(r16)
#define Dpx0 24(r16)
#define Dpx1 28(r16)
#define Dpx2 32(r16)
#define Dpy0 36(r16)
#define Dpy1 40(r16)
#define Dpy2 44(r16)
#define Dpz0 48(r16)
#define Dpz1 52(r16)
#define Dpz2 56(r16)
#define Qpxx0 60(r16)
#define Qpxx1 64(r16)
#define Qpxx2 68(r16)
#define Qpyy0 72(r16)
#define Qpyy1 76(r16)
#define Qpyy2 80(r16)
#define Qpzz0 84(r16)
#define Qpzz1 88(r16)
#define Qpzz2 92(r16)
#define Qpxy0 96(r16)
#define Qpxy1 100(r16)
#define Qpxy2 104(r16)
#define Qpxz0 108(r16)
#define Qpxz1 112(r16)
#define Qpxz2 116(r16)
#define Qpyz0 120(r16)
#define Qpyz1 124(r16)
#define Qpyz2 128(r16)

#define wx0 32(sp)
#define wx1 36(sp)
#define wx2 40(sp)

	.text
	.align 32
_InteractCell::
.a2 = 0
.f2 = 72
	addu -(.a2+.f2), sp, sp
	pfmul.ss f0, f0, f0	// these six lines result in 6x speedup
	pfmul.ss f0, f0, f0	// if there is junk in pipes.  Why?
	pfmul.ss f0, f0, f0
	pfadd.ss f0, f0, f0
	pfadd.ss f0, f0, f0
	pfadd.ss f0, f0, f0
					fld.l 0(r17),x
					fld.l 4(r17),y
.dual
	r2apt.ss f0,f0,f0;		fld.l 8(r17),z
	r2apt.ss f0,f0,f0;		fld.l 0(r16),t0
	r2apt.ss f0,f0,f0;		fld.l 4(r16),t1
	i2apt.ss f0,f0,f0;		fld.l 8(r16),t2
	pfsub.ss x, t0, f0;		fst.d f2, 8(sp)
.enddual
	pfsub.ss y, t1, f0;		fst.d f4, 16(sp)
	pfsub.ss z, t2, f0;		fst.d f6, 24(sp)
	pfsub.ss f0, f0, x		
	pfsub.ss f0, f0, y		
	pfsub.ss f0, f0, z		
	pfmul.ss x, x, f0
	pfmul.ss y, y, f0
	pfmul.ss z, z, f0
	pfmul.ss f0, f0, xx
	pfmul.ss f0, f0, yy
	pfmul.ss f0, f0, zz
	fadd.ss xx,yy,t0
		
	
	fadd.ss zz, t0, f27

	// inline recipsqrtf
        orh 0x3f00, r0, r29
        frsqr.ss f27, t1
        ixfr r29, t0
        fmul.ss t1, t1, f29
        orh 0x3fc0, r0, r29
        ixfr r29, f28
        fmul.ss t0, f27, t0
        
	
	fmul.ss f29, t0, f29
        
	
	fsub.ss f28, f29, f29
        
	
	fmul.ss t1, f29, t1
					orh h%_CellCnt+4,r0,r29
					or l%_CellCnt+4,r29,r29
	fmul.ss t1, t1, f29
					ld.l 0(r29),r30
					addu 1,r30,r30
	fmul.ss t0, f29, f29
					st.l r30,0(r29)
					ld.l 156(r17),r30
	fsub.ss f28, f29, f29
					addu 1,r30,r30		// 1 nterm
					st.l r30,156(r17)
	fmul.ss t1, f29, cf0
					orh h%0x40400000,r0,r29
					ixfr r29,t0
	fmul.ss cf0,cf0,f22    // d2inv
					orh h%0x40a00000,r0,r29
					ixfr r29,t1
	fmul.ss f22,cf0,cf1
					orh h%0x40e00000,r0,r29
					ixfr r29,t2
	fmul.ss f22,t0,t0


	fmul.ss f22,t1,t1


	fmul.ss f22,t2,t2


	fmul.ss t0,cf1,cf3


	fmul.ss t1,cf3,cf5    // cf5
					orh h%0x3f000000,r0,r20
					ixfr r20,f20
	fmul.ss t2,cf5,t2    // cf7

.dual
	pfmul.ss cf3,f20,f0;		fld.l Dpx0, t0
	pfmul.ss cf5,f20,f0;		fld.l Dpx1, t1
	pfmul.ss t2, f20,f0;		fld.l Dpx2, t2
	// Ddot(XdotD, x, y, z, Dpole(cp).x, Dpole(cp).y, Dpole(cp).z)
	pfmul.ss t0, x, cf3_2;		fld.l Dpy0, t0
	pfmul.ss t1, x, cf5_2;		fld.l Dpy1, t1
	pfmul.ss t2, x, cf7_2;		fld.l Dpy2, t2
	m12tpm.ss t0, y, f0;		fld.l Dpz0, t0
	m12tpm.ss t1, y, f0;		fld.l Dpz1, t1
	m12tpm.ss t2, y, f0;		fld.l Dpz2, t2
	m12apm.ss t0, z, f0;		fld.l Qpxx0, t0
	m12apm.ss t1, z, f0;		fld.l Qpxx1, t1
	m12apm.ss t2, z, f0;		fld.l Qpxx2, t2
	// Ddot(QdotX.x, x, y, z, Qpole(cp).xx, Qpole(cp).xy, Qpole(cp).xz)
	m12apm.ss t0, x, f0;		fld.l Qpxy0, t0
	m12apm.ss t1, x, f0;		fld.l Qpxy1, t1
	m12apm.ss t2, x, f0;		fld.l Qpxy2, t2
	m12tpm.ss t0, y, Xd0;		fld.l Qpxz0, t0
	m12tpm.ss t1, y, Xd1;		fld.l Qpxz1, t1
	m12tpm.ss t2, y, Xd2;		fld.l Qpxz2, t2
	m12apm.ss t0, z, f0;		fld.l Qpxy0, t0
	m12apm.ss t1, z, f0;		fld.l Qpxy1, t1
	m12apm.ss t2, z, f0;		fld.l Qpxy2, t2
	// Ddot(QdotX.y, x, y, z, Qpole(cp).xy, Qpole(cp).yy, Qpole(cp).yz)
	m12apm.ss t0, x, f0;		fld.l Qpyy0, t0
	m12apm.ss t1, x, f0;		fld.l Qpyy1, t1
	m12apm.ss t2, x, f0;		fld.l Qpyy2, t2
	m12tpm.ss t0, y, Qx0;		fld.l Qpyz0, t0
	m12tpm.ss t1, y, Qx1;		fld.l Qpyz1, t1
	m12tpm.ss t2, y, Qx2;		fld.l Qpyz2, t2
	m12apm.ss t0, z, f0;		fld.l Qpxz0, t0
	m12apm.ss t1, z, f0;		fld.l Qpxz1, t1
	m12apm.ss t2, z, f0;		fld.l Qpxz2, t2
	// Ddot(QdotX.z, x, y, z, Qpole(cp).xz, Qpole(cp).yz, Qpole(cp).zz)
	m12apm.ss t0, x, f0;		fld.l Qpyz0, t0
	m12apm.ss t1, x, f0;		fld.l Qpyz1, t1
	m12apm.ss t2, x, f0;		fld.l Qpyz2, t2
	m12tpm.ss t0, y, Qy0;		fld.l Qpzz0, t0
	m12tpm.ss t1, y, Qy1;		fld.l Qpzz1, t1
	m12tpm.ss t2, y, Qy2;		fld.l Qpzz2, t2
	m12apm.ss t0, z, f0;		fld.l St0, st0
	m12apm.ss t1, z, f0;		fld.l St1, st1
	m12apm.ss t2, z, f0;		fld.l St2, st2
	// Ddot(XdotQdotX, x, y, z, QdotX.x, QdotX.y, QdotX.z)
	m12apm.ss Qx0, x, f0;		fld.l 28(r17), t0	// load psi
	m12apm.ss Qx1, x, f0;		fld.l 32(r17), t1
.enddual
	m12apm.ss Qx2, x, f0;		fld.l 36(r17), t2
	m12tpm.ss Qy0, y, Qz0;		nop
	m12tpm.ss Qy1, y, Qz1
	m12tpm.ss Qy2, y, Qz2
	m12apm.ss Qz0, z, f0
	m12apm.ss Qz1, z, f0
	m12apm.ss Qz2, z, f0
	// Ddot(tmppsi, cf0, cf1, 0.5F*cf3, Strength(cp), XdotD, XdotQdotX)
	m12apm.ss st0, cf0, f0
	m12apm.ss st1, cf0, f0
	m12apm.ss st2, cf0, f0
	pfadd.ss f0, t0, XQ0
	pfadd.ss f0, t1, XQ1
	pfadd.ss f0, t2, XQ2
	m12apm.ss Xd0, cf1, f0
	m12apm.ss Xd1, cf1, f0
	m12apm.ss Xd2, cf1, f0
	m12apm.ss XQ0, cf3_2, f0
	m12apm.ss XQ1, cf3_2, f0
	m12apm.ss XQ2, cf3_2, f0
	// Ddot(tmpa, cf1, cf3, 0.5F*cf5, Strength(cp), XdotD, XdotQdotX)
	m12apm.ss st0, cf1, f0
	m12apm.ss st1, cf1, f0
	m12apm.ss st2, cf1, f0
	m12tpm.ss Xd0, cf3, t0
	m12tpm.ss Xd1, cf3, t1
	m12tpm.ss Xd2, cf3, t2
	m12apm.ss XQ0, cf5_2, f0
	m12apm.ss XQ1, cf5_2, f0
	m12apm.ss XQ2, cf5_2, f0
	// Ddot(tmpb, cf3, cf5, 0.5F*cf7, Strength(cp), XdotD, XdotQdotX)
	m12apm.ss st0, cf3, f0
	m12apm.ss st1, cf3, f0
	m12apm.ss st2, cf3, f0
	m12tpm.ss Xd0, cf5, ta0
.dual
	m12tpm.ss Xd1, cf5, ta1
	m12tpm.ss Xd2, cf5, ta2
	m12apm.ss XQ0, cf7_2, f0;	fst.l t0, 28(r17)	// store psi
	m12apm.ss XQ1, cf7_2, f0;	fst.l t1, 32(r17)
	m12apm.ss XQ2, cf7_2, f0;	fst.l t2, 36(r17)

	// Vel(me)[0] += -y*tmpa[2]+cf1*Dpole(cp).y[2]+cf3*QdotX.y[2];
	// Vel(me)[1] += -z*tmpa[0]+cf1*Dpole(cp).z[0]+cf3*QdotX.z[0];
	// Vel(me)[2] += -x*tmpa[1]+cf1*Dpole(cp).x[1]+cf3*QdotX.x[1];

	// Vel(me)[0] -= -z*tmpa[1]+cf1*Dpole(cp).z[1]+cf3*QdotX.z[1];
	// Vel(me)[1] -= -x*tmpa[2]+cf1*Dpole(cp).x[2]+cf3*QdotX.x[2];
	// Vel(me)[2] -= -y*tmpa[0]+cf1*Dpole(cp).y[0]+cf3*QdotX.y[0];

	m12apm.ss ta2, y, f0;		fld.l 40(r17), t0	// load vel
	m12apm.ss ta0, z, f0;		fld.l 44(r17), t1
	m12apm.ss ta1, x, f0;		fld.l 48(r17), t2
	pfadd.ss f0, t0, tb0;		fld.l Dpy2, t0
	pfadd.ss f0, t1, tb1;		fld.l Dpz0, t1
	pfadd.ss f0, t2, tb2;		fld.l Dpx1, t2
	m12asm.ss t0, cf1, f0;		nop
	m12asm.ss t1, cf1, f0;		nop
	m12asm.ss t2, cf1, f0;		nop
	m12apm.ss Qy2, cf3, f0;		nop
	m12apm.ss Qz0, cf3, f0;		nop
	m12apm.ss Qx1, cf3, f0;		nop
	m12apm.ss ta1, z, f0;		fld.l Dpz1, t0
	m12apm.ss ta2, x, f0;		fld.l Dpx2, t1
	m12apm.ss ta0, y, f0;		fld.l Dpy0, t2
	m12apm.ss t0, cf1, f0;		nop
	m12apm.ss t1, cf1, f0;		nop
	m12apm.ss t2, cf1, f0;		nop
	m12asm.ss Qz1, cf3, f0;		fld.l Dpx0, t0
	m12asm.ss Qx2, cf3, f0;		fld.l Dpx1, t1
	m12asm.ss Qy0, cf3, f0;		fld.l Dpx2, t2
	// Ddot2(wrk.x, cf3, cf5, Dpole(cp).x, QdotX.x)
	m12asm.ss t0, cf3, f0;		fld.l Dpy0, t0
	m12asm.ss t1, cf3, f0;		fld.l Dpy1, t1
	m12asm.ss t2, cf3, f0;		fld.l Dpy2, t2
	m12tpm.ss Qx0, cf5, s0;		fst.l s0, 40(r17)	// store vel
	m12tpm.ss Qx1, cf5, s1;		fst.l s1, 44(r17)
	m12tpm.ss Qx2, cf5, s2;		fst.l s2, 48(r17)
	// Ddot2(wrk.y, cf3, cf5, Dpole(cp).y, QdotX.y)
	m12apm.ss t0, cf3, f0;		fld.l Dpz0, t0
	m12apm.ss t1, cf3, f0;		fld.l Dpz1, t1
	m12apm.ss t2, cf3, f0;		fld.l Dpz2, t2
	m12tpm.ss Qy0, cf5, s0;		fst.l s0, wx0
	m12tpm.ss Qy1, cf5, s1;		fst.l s1, wx1
	m12tpm.ss Qy2, cf5, s2;		fst.l s2, wx2
	// Ddot(wrk.z, cf3, cf5, Dpole(cp).z, QdotX.z)
	m12apm.ss t0, cf3, f0;		nop
	m12apm.ss t1, cf3, f0;		nop
	m12apm.ss t2, cf3, f0;		nop
	m12tpm.ss Qz0, cf5, wy0;	nop
	m12tpm.ss Qz1, cf5, wy1;	nop
	m12tpm.ss Qz2, cf5, wy2;	nop
	m12apm.ss x,x,f0;  		nop
	m12apm.ss y,y,f0;		nop
	m12apm.ss z,z,f0;		nop
	pfmul.ss x,y,xx;		orh h%0xc0000000,r0,r20
	pfmul.ss x,z,yy;		ixfr r20,t2
	pfmul.ss y,z,zz;		nop
	pfmul.ss x,t2,xy;		fld.l 76(r17),t0
	pfmul.ss y,t2,xz;		fld.l 56(r17),t1
	pfmul.ss z,t2,yz;		fld.l 72(r17),t2
	// Ddot(dxxpsi, x*x, -2.0F*x, cf3, tmpb, wrk.x, Qpole(cp).xx, -tmpa)
				// load gradvel, store wrk.z
	pfsub.ss t0, ta1, wz0; 		nop
	pfsub.ss t1, ta2, wz1;          nop
	pfsub.ss t2, ta0, wz2;          nop
	pfmul.ss tb1, xx, m2x;		fld.l wx1, t0
	pfmul.ss tb2, yy, m2y;		nop
	pfmul.ss tb0, zz, m2z;		nop
	m12apm.ss t0, m2x, f0;		fld.l Qpxx1, t0
	m12apm.ss wy2, m2y, f0;		fld.l Qpyy2, t1
	m12apm.ss wz0, m2z, f0;		fld.l Qpzz0, t2
	m12apm.ss t0, cf3, f0;		nop
	m12apm.ss t1, cf3, f0;		nop
	m12apm.ss t2, cf3, f0;		nop
	m12apm.ss tb2, xx, f0;		fld.l 64(r17), t0
	m12apm.ss tb0, yy, f0;		fld.l 80(r17), t1
	m12apm.ss tb1, zz, f0;		fld.l 60(r17), t2
	pfadd.ss t0, ta2, xxp1;		fld.l wx2, t0
	pfadd.ss t1, ta0, yyp2;		nop
	pfadd.ss t2, ta1, zzp0;		nop
	m12asm.ss t0, m2x, f0;		fld.l Qpxx2, t0
	m12asm.ss wy0, m2y, f0;		fld.l Qpyy0, t1
	m12asm.ss wz1, m2z, f0;		fld.l Qpzz1, t2
	m12asm.ss t0, cf3, f0;		nop
	m12asm.ss t1, cf3, f0;		nop
	m12asm.ss t2, cf3, f0;		nop
	// Ddot(dxypsi, x*y, -x, -y, cf3, tmpb, wrk.y, wrk.x, Qpole(cp).xy)
	m12asm.ss tb0, xy, f0;		nop
	m12asm.ss tb1, xy, f0;		nop
	m12asm.ss tb2, xy, f0;		nop
	m12tpm.ss wy0, x, xxp2;		fld.l wx0, t0
	m12tpm.ss wy1, x, yyp0;		fld.l wx1, t1
	m12tpm.ss wy2, x, zzp1;		fld.l wx2, t2
	m12asm.ss t0, y, f0;		fld.l Qpxy0, t0
	m12asm.ss t1, y, f0;		fld.l Qpxy1, t1
	m12asm.ss t2, y, f0;		fld.l Qpxy2, t2
	m12asm.ss t0, cf3, f0;		nop
	m12asm.ss t1, cf3, f0;		nop
	m12asm.ss t2, cf3, f0;		nop
	// Ddot(dxzpsi, x*z, -x, -z, cf3, tmpb, wrk.z, wrk.x, Qpole(cp).xz)
	m12apm.ss tb0, xz, f0;		nop
	m12apm.ss tb1, xz, f0;		nop
	m12apm.ss tb2, xz, f0;		nop
	m12tpm.ss wz0, x, xyp0;		fld.l wx0, t0
	m12tpm.ss wz1, x, xyp1;		fld.l wx1, t1
	m12tpm.ss wz2, x, xyp2;		fld.l wx2, t2
	m12asm.ss t0, z, f0;		fld.l Qpxz0, t0
	m12asm.ss t1, z, f0;		fld.l Qpxz1, t1
	m12asm.ss t2, z, f0;		fld.l Qpxz2, t2	
	m12asm.ss t0, cf3, f0;		nop
	m12asm.ss t1, cf3, f0;		nop
	m12asm.ss t2, cf3, f0;		nop
	// Ddot(dyzpsi, y*z, -y, -z, cf3, tmpb, wrk.z, wrk.y, Qpole(cp).yz)
	m12apm.ss tb0, yz, f0;		nop
	m12apm.ss tb1, yz, f0;		nop
	m12apm.ss tb2, yz, f0;		nop
	m12tpm.ss wz0, y, xzp0;		nop
	m12tpm.ss wz1, y, xzp1;		nop
	m12tpm.ss wz2, y, xzp2;		nop
	m12asm.ss wy0, z, f0;		fld.l Qpyz0, t0
	m12asm.ss wy1, z, f0;		fld.l Qpyz1, t1
	m12asm.ss wy2, z, f0;		fld.l Qpyz2, t2
	m12asm.ss t0, cf3, f0;		fld.l 52(r17),t0
	m12asm.ss t1, cf3, f0;		fld.l 68(r17),t1
	m12asm.ss t2, cf3, f0;		fld.l 84(r17),t2
	m12apm.ss f0, f0, f0;		nop
	m12apm.ss f0, f0, f0;		nop
	m12apm.ss f0, f0, f0;		nop
	pfadd.ss xyp2, t0, yzp0;	fld.d 8(sp), f2
	pfadd.ss yzp0, t1, yzp1;	fld.d 16(sp), f4
	pfadd.ss xzp1, t2, yzp2;	fld.d 24(sp), f6
	pfsub.ss yyp2,yzp1,x;		nop
	pfadd.ss yzp2,zzp1,y;		nop
	pfadd.ss xzp0,xxp2,z;		nop
	pfsub.ss x, xzp1, t0;		fst.l t0,56(r17)
	pfsub.ss y, xyp2, t1;		fst.l t1,60(r17)
	pfsub.ss z, yzp0, t2;		fst.l t2,64(r17)
	pfsub.ss zzp0,xzp2,t0;		fst.l t0,52(r17)
	pfsub.ss xxp1,xyp0,t1;		fst.l t1,68(r17)
	pfadd.ss xyp1,yyp0,t2;		fst.l t2,84(r17)
	pfadd.ss f0, f0, t0;		fst.l t0,72(r17)
.enddual
	pfadd.ss f0, f0, t1;		fst.l t1,76(r17)
	pfadd.ss f0, f0, t2;		fst.l t2,80(r17)
	bri r1
	addu (.a2+.f2), sp, sp
.comm _CellCnt,40
