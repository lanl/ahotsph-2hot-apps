	.globl	_clear_tregs
	.align	32
_clear_tregs:
	r2apt.ss f0,f0,f0
	r2apt.ss f0,f0,f0
	r2apt.ss f0,f0,f0
	i2apt.ss f0,f0,f0
	bri r1
	nop

// void Ddot(float *z, float r0, float r1, float r2, 
//	  const float *y0, const float *y1, const float *y2);

	.globl	_Ddot_asm
	.align	32
_Ddot_asm:
					fld.l 0(r17), f16
					fld.l 4(r17), f17
	d.pfmul.ss f0, f0, f0	// dummy instruction to initiate dual mode
					fld.l 8(r17), f18
	d.pfmul.ss f16, f8, f0  // enter dual mode
					fld.l 0(r18), f16
	d.pfmul.ss f17, f8, f0
					fld.l 4(r18), f17
	d.pfmul.ss f18, f8, f0
					fld.l 8(r18), f18
	d.m12tpm.ss f16, f9, f0
					fld.l 0(r19), f16
	m12tpm.ss f17, f9, f0
					fld.l 4(r19), f17
	m12tpm.ss f18, f9, f0
					fld.l 8(r19), f18   // exit dual mode
	m12apm.ss f16, f10, f0
	m12apm.ss f17, f10, f0
	m12apm.ss f18, f10, f0
	m12apm.ss f0, f0, f0
	d.m12apm.ss f0, f0, f0
	d.m12apm.ss f0, f0, f0
	d.pfadd.ss f0, f0, f16  // enter dual mode
					fst.l f16, r0(r16)
	pfadd.ss f0, f0, f17
					fst.l f17, 4(r16)
	pfadd.ss f0, f0, f18
					fst.l f18, 8(r16)   // exit dual mode
	bri r1
	nop

	.globl	_Ddot_asm_nodual
	.align	32
_Ddot_asm_nodual:	
	fld.l 0(r17), f16
	fld.l 4(r17), f17
	fld.l 8(r17), f18
	pfmul.ss f16, f8, f0
	fld.l 0(r18), f16
	pfmul.ss f17, f8, f0
	fld.l 4(r18), f17
	pfmul.ss f18, f8, f0
	fld.l 8(r18), f18
	m12tpm.ss f16, f9, f0
	fld.l 0(r19), f16
	m12tpm.ss f17, f9, f0
	fld.l 4(r19), f17
	m12tpm.ss f18, f9, f0
	fld.l 8(r19), f18
	m12apm.ss f16, f10, f0
	m12apm.ss f17, f10, f0
	m12apm.ss f18, f10, f0
	m12apm.ss f0, f0, f0
	m12apm.ss f0, f0, f0
	m12apm.ss f0, f0, f0
	pfadd.ss f0, f0, f16
	fst.l f16, r0(r16)
	pfadd.ss f0, f0, f17
	fst.l f17, 4(r16)
	pfadd.ss f0, f0, f18
	fst.l f18, 8(r16)
	bri r1
	nop
