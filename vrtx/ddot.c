
void Ddot3(float *z, float r0, float r1, float r2, const float *y);

main(int argc, char *argv[])
{
    float z[3];
    float y[9] = {1., 2., 3., 4., 5., 6., 7., 8., 9.};

    Ddot3(z, 1., 2., 3., y);

    printf("%f %f %f\n", z[0], z[1], z[2]);
      
    exit(0);
}

	.globl	_Ddot3
	.align	32
_Ddot3:
.dual
	pfadd.ss f0, f0, f0
					fld.l 0(r17), f16
	pfadd.ss f0, f0, f0
					fld.l 4(r17), f17
	pfadd.ss f0, f0, f0
					fld.l 8(r17), f18
	pfmul.ss f16, f8, f0
					fld.l 12(r17), f19
	pfmul.ss f17, f9, f0
					fld.l 16(r17), f16
	pfmul.ss f18, f10, f0
					fld.l 20(r17), f17
	m12apm.ss f19, f8, f0
					fld.l 24(r17), f18
	m12apm.ss f16, f9, f0
					fld.l 28(r17), f19
	m12apm.ss f17, f10, f0
					fld.l 32(r17), f16
	m12apm.ss f18, f8, f0
					nop
	m12apm.ss f19, f9, f0
					nop
	m12apm.ss f16, f10, f0
					nop
	pfadd.ss f0, f0, f16
					fst.l f16, r0(r16)
	pfadd.ss f0, f0, f17
					fst.l f17, 4(r16)
.enddual	
	pfadd.ss f0, f0, f18
	bri r1 
	fst.l f18, 8(r16)
