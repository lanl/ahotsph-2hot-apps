/* rfftpack.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $$$Received from netlib in response to */
/* $$$Subject: send rffti rfftf rfftb from fftpack */
/*$$$Caveat receptor.  (Jack) dongarra@cs.utk.edu, (Eric Grosse) research!ehg
*/
/* $$$Careful! Anything free comes with no guarantee. */
/* $$$* */
/* $$$*** from netlib, Mon Dec 10 17:35:52 EST 1990 *** */
/* Subroutine */ int rfftb_(n, r, wsave)
integer *n;
real *r, *wsave;
{
    extern /* Subroutine */ int rfftb1_();

    /* Parameter adjustments */
    --wsave;
    --r;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    rfftb1_(n, &r[1], &wsave[1], &wsave[*n + 1], &wsave[(*n << 1) + 1]);
    return 0;
} /* rfftb_ */

/* Subroutine */ int rfftb1_(n, c, ch, wa, ifac)
integer *n;
real *c, *ch, *wa;
integer *ifac;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int radb2_(), radb3_(), radb4_(), radb5_();
    static integer i;
    extern /* Subroutine */ int radbg_();
    static integer k1, l1, l2, na, nf, ip, iw, ix2, ix3, ix4, ido, idl1;

    /* Parameter adjustments */
    --ifac;
    --wa;
    --ch;
    --c;

    /* Function Body */
    nf = ifac[2];
    na = 0;
    l1 = 1;
    iw = 1;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idl1 = ido * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	if (na != 0) {
	    goto L101;
	}
	radb4_(&ido, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	radb4_(&ido, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;
L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	radb2_(&ido, &l1, &c[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	radb2_(&ido, &l1, &ch[1], &c[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;
L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + ido;
	if (na != 0) {
	    goto L107;
	}
	radb3_(&ido, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	radb3_(&ido, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;
L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	ix4 = ix3 + ido;
	if (na != 0) {
	    goto L110;
	}
	radb5_(&ido, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]
		);
	goto L111;
L110:
	radb5_(&ido, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]
		);
L111:
	na = 1 - na;
	goto L115;
L112:
	if (na != 0) {
	    goto L113;
	}
	radbg_(&ido, &ip, &l1, &idl1, &c[1], &c[1], &c[1], &ch[1], &ch[1], &
		wa[iw]);
	goto L114;
L113:
	radbg_(&ido, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c[1], &c[1], &
		wa[iw]);
L114:
	if (ido == 1) {
	    na = 1 - na;
	}
L115:
	l1 = l2;
	iw += (ip - 1) * ido;
/* L116: */
    }
    if (na == 0) {
	return 0;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	c[i] = ch[i];
/* L117: */
    }
    return 0;
} /* rfftb1_ */

/* Subroutine */ int rfftf_(n, r, wsave)
integer *n;
real *r, *wsave;
{
    extern /* Subroutine */ int rfftf1_();

    /* Parameter adjustments */
    --wsave;
    --r;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    rfftf1_(n, &r[1], &wsave[1], &wsave[*n + 1], &wsave[(*n << 1) + 1]);
    return 0;
} /* rfftf_ */

/* Subroutine */ int rfftf1_(n, c, ch, wa, ifac)
integer *n;
real *c, *ch, *wa;
integer *ifac;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int radf2_(), radf3_(), radf4_(), radf5_();
    static integer i;
    extern /* Subroutine */ int radfg_();
    static integer k1, l1, l2, na, kh, nf, ip, iw, ix2, ix3, ix4, ido, idl1;

    /* Parameter adjustments */
    --ifac;
    --wa;
    --ch;
    --c;

    /* Function Body */
    nf = ifac[2];
    na = 1;
    l2 = *n;
    iw = *n;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	kh = nf - k1;
	ip = ifac[kh + 3];
	l1 = l2 / ip;
	ido = *n / l2;
	idl1 = ido * l1;
	iw -= (ip - 1) * ido;
	na = 1 - na;
	if (ip != 4) {
	    goto L102;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	if (na != 0) {
	    goto L101;
	}
	radf4_(&ido, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L110;
L101:
	radf4_(&ido, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L110;
L102:
	if (ip != 2) {
	    goto L104;
	}
	if (na != 0) {
	    goto L103;
	}
	radf2_(&ido, &l1, &c[1], &ch[1], &wa[iw]);
	goto L110;
L103:
	radf2_(&ido, &l1, &ch[1], &c[1], &wa[iw]);
	goto L110;
L104:
	if (ip != 3) {
	    goto L106;
	}
	ix2 = iw + ido;
	if (na != 0) {
	    goto L105;
	}
	radf3_(&ido, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L110;
L105:
	radf3_(&ido, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2]);
	goto L110;
L106:
	if (ip != 5) {
	    goto L108;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	ix4 = ix3 + ido;
	if (na != 0) {
	    goto L107;
	}
	radf5_(&ido, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]
		);
	goto L110;
L107:
	radf5_(&ido, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]
		);
	goto L110;
L108:
	if (ido == 1) {
	    na = 1 - na;
	}
	if (na != 0) {
	    goto L109;
	}
	radfg_(&ido, &ip, &l1, &idl1, &c[1], &c[1], &c[1], &ch[1], &ch[1], &
		wa[iw]);
	na = 1;
	goto L110;
L109:
	radfg_(&ido, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c[1], &c[1], &
		wa[iw]);
	na = 0;
L110:
	l2 = l1;
/* L111: */
    }
    if (na == 1) {
	return 0;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	c[i] = ch[i];
/* L112: */
    }
    return 0;
} /* rfftf1_ */

/* Subroutine */ int rffti_(n, wsave)
integer *n;
real *wsave;
{
    extern /* Subroutine */ int rffti1_();

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    rffti1_(n, &wsave[*n + 1], &wsave[(*n << 1) + 1]);
    return 0;
} /* rffti_ */

/* Subroutine */ int rffti1_(n, wa, ifac)
integer *n;
real *wa;
integer *ifac;
{
    /* Initialized data */

    static integer ntryh[4] = { 4,2,3,5 };

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double cos(), sin();

    /* Local variables */
    static real argh;
    static integer ntry, i, j;
    static real argld;
    static integer k1, l1, l2, ib;
    static real fi;
    static integer ld, ii, nf, ip, nl, is, nq, nr;
    static real arg;
    static integer ido, ipm;
    static real tpi;
    static integer nfm1;

    /* Parameter adjustments */
    --ifac;
    --wa;

    /* Function Body */
    nl = *n;
    nf = 0;
    j = 0;
L101:
    ++j;
    if (j - 4 <= 0) {
	goto L102;
    } else {
	goto L103;
    }
L102:
    ntry = ntryh[j - 1];
    goto L104;
L103:
    ntry += 2;
L104:
    nq = nl / ntry;
    nr = nl - ntry * nq;
    if (nr != 0) {
	goto L101;
    } else {
	goto L105;
    }
L105:
    ++nf;
    ifac[nf + 2] = ntry;
    nl = nq;
    if (ntry != 2) {
	goto L107;
    }
    if (nf == 1) {
	goto L107;
    }
    i__1 = nf;
    for (i = 2; i <= i__1; ++i) {
	ib = nf - i + 2;
	ifac[ib + 2] = ifac[ib + 1];
/* L106: */
    }
    ifac[3] = 2;
L107:
    if (nl != 1) {
	goto L104;
    }
    ifac[1] = *n;
    ifac[2] = nf;
    tpi = (float)6.28318530717959;
    argh = tpi / (real) (*n);
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;
    if (nfm1 == 0) {
	return 0;
    }
    i__1 = nfm1;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	ipm = ip - 1;
	i__2 = ipm;
	for (j = 1; j <= i__2; ++j) {
	    ld += l1;
	    i = is;
	    argld = (real) ld * argh;
	    fi = (float)0.;
	    i__3 = ido;
	    for (ii = 3; ii <= i__3; ii += 2) {
		i += 2;
		fi += (float)1.;
		arg = fi * argld;
		wa[i - 1] = cos(arg);
		wa[i] = sin(arg);
/* L108: */
	    }
	    is += ido;
/* L109: */
	}
	l1 = l2;
/* L110: */
    }
    return 0;
} /* rffti1_ */

/* Subroutine */ int radb2_(ido, l1, cc, ch, wa1)
integer *ido, *l1;
real *cc, *ch, *wa1;
{
    /* System generated locals */
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i, k, ic;
    static real ti2, tr2;
    static integer idp2;

    /* Parameter adjustments */
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 3 + 1;
    cc -= cc_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] + 
		cc[*ido + ((k << 1) + 2) * cc_dim1];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 
		+ 1] - cc[*ido + ((k << 1) + 2) * cc_dim1];
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = cc[i - 1 + ((k << 1) + 1) * 
		    cc_dim1] + cc[ic - 1 + ((k << 1) + 2) * cc_dim1];
	    tr2 = cc[i - 1 + ((k << 1) + 1) * cc_dim1] - cc[ic - 1 + ((k << 1)
		     + 2) * cc_dim1];
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + ((k << 1) + 1) * cc_dim1]
		     - cc[ic + ((k << 1) + 2) * cc_dim1];
	    ti2 = cc[i + ((k << 1) + 1) * cc_dim1] + cc[ic + ((k << 1) + 2) * 
		    cc_dim1];
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 2] * tr2 - 
		    wa1[i - 1] * ti2;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 2] * ti2 + wa1[i 
		    - 1] * tr2;
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[*ido + (k + ch_dim2) * ch_dim1] = cc[*ido + ((k << 1) + 1) * 
		cc_dim1] + cc[*ido + ((k << 1) + 1) * cc_dim1];
	ch[*ido + (k + (ch_dim2 << 1)) * ch_dim1] = -(doublereal)(cc[((k << 1)
		 + 2) * cc_dim1 + 1] + cc[((k << 1) + 2) * cc_dim1 + 1]);
/* L106: */
    }
L107:
    return 0;
} /* radb2_ */

/* Subroutine */ int radb3_(ido, l1, cc, ch, wa1, wa2)
integer *ido, *l1;
real *cc, *ch, *wa1, *wa2;
{
    /* Initialized data */

    static real taur = (float)-.5;
    static real taui = (float).866025403784439;

    /* System generated locals */
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i, k, ic;
    static real ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
    static integer idp2;

    /* Parameter adjustments */
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = (cc_dim1 << 2) + 1;
    cc -= cc_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	tr2 = cc[*ido + (k * 3 + 2) * cc_dim1] + cc[*ido + (k * 3 + 2) * 
		cc_dim1];
	cr2 = cc[(k * 3 + 1) * cc_dim1 + 1] + taur * tr2;
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 3 + 1) * cc_dim1 + 1] + tr2;
	ci3 = taui * (cc[(k * 3 + 3) * cc_dim1 + 1] + cc[(k * 3 + 3) * 
		cc_dim1 + 1]);
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr2 + ci3;
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    tr2 = cc[i - 1 + (k * 3 + 3) * cc_dim1] + cc[ic - 1 + (k * 3 + 2) 
		    * cc_dim1];
	    cr2 = cc[i - 1 + (k * 3 + 1) * cc_dim1] + taur * tr2;
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = cc[i - 1 + (k * 3 + 1) * 
		    cc_dim1] + tr2;
	    ti2 = cc[i + (k * 3 + 3) * cc_dim1] - cc[ic + (k * 3 + 2) * 
		    cc_dim1];
	    ci2 = cc[i + (k * 3 + 1) * cc_dim1] + taur * ti2;
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * 3 + 1) * cc_dim1] + 
		    ti2;
	    cr3 = taui * (cc[i - 1 + (k * 3 + 3) * cc_dim1] - cc[ic - 1 + (k *
		     3 + 2) * cc_dim1]);
	    ci3 = taui * (cc[i + (k * 3 + 3) * cc_dim1] + cc[ic + (k * 3 + 2) 
		    * cc_dim1]);
	    dr2 = cr2 - ci3;
	    dr3 = cr2 + ci3;
	    di2 = ci2 + cr3;
	    di3 = ci2 - cr3;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 2] * dr2 - 
		    wa1[i - 1] * di2;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 2] * di2 + wa1[i 
		    - 1] * dr2;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 2] * dr3 - wa2[
		    i - 1] * di3;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 2] * di3 + wa2[i - 
		    1] * dr3;
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* radb3_ */

/* Subroutine */ int radb4_(ido, l1, cc, ch, wa1, wa2, wa3)
integer *ido, *l1;
real *cc, *ch, *wa1, *wa2, *wa3;
{
    /* Initialized data */

    static real sqrt2 = (float)1.414213562373095;

    /* System generated locals */
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i, k, ic;
    static real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, 
	    tr3, tr4;
    static integer idp2;

    /* Parameter adjustments */
    --wa3;
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 5 + 1;
    cc -= cc_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	tr1 = cc[((k << 2) + 1) * cc_dim1 + 1] - cc[*ido + ((k << 2) + 4) * 
		cc_dim1];
	tr2 = cc[((k << 2) + 1) * cc_dim1 + 1] + cc[*ido + ((k << 2) + 4) * 
		cc_dim1];
	tr3 = cc[*ido + ((k << 2) + 2) * cc_dim1] + cc[*ido + ((k << 2) + 2) *
		 cc_dim1];
	tr4 = cc[((k << 2) + 3) * cc_dim1 + 1] + cc[((k << 2) + 3) * cc_dim1 
		+ 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = tr2 + tr3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = tr1 - tr4;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = tr2 - tr3;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = tr1 + tr4;
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    ti1 = cc[i + ((k << 2) + 1) * cc_dim1] + cc[ic + ((k << 2) + 4) * 
		    cc_dim1];
	    ti2 = cc[i + ((k << 2) + 1) * cc_dim1] - cc[ic + ((k << 2) + 4) * 
		    cc_dim1];
	    ti3 = cc[i + ((k << 2) + 3) * cc_dim1] - cc[ic + ((k << 2) + 2) * 
		    cc_dim1];
	    tr4 = cc[i + ((k << 2) + 3) * cc_dim1] + cc[ic + ((k << 2) + 2) * 
		    cc_dim1];
	    tr1 = cc[i - 1 + ((k << 2) + 1) * cc_dim1] - cc[ic - 1 + ((k << 2)
		     + 4) * cc_dim1];
	    tr2 = cc[i - 1 + ((k << 2) + 1) * cc_dim1] + cc[ic - 1 + ((k << 2)
		     + 4) * cc_dim1];
	    ti4 = cc[i - 1 + ((k << 2) + 3) * cc_dim1] - cc[ic - 1 + ((k << 2)
		     + 2) * cc_dim1];
	    tr3 = cc[i - 1 + ((k << 2) + 3) * cc_dim1] + cc[ic - 1 + ((k << 2)
		     + 2) * cc_dim1];
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = tr2 + tr3;
	    cr3 = tr2 - tr3;
	    ch[i + (k + ch_dim2) * ch_dim1] = ti2 + ti3;
	    ci3 = ti2 - ti3;
	    cr2 = tr1 - tr4;
	    cr4 = tr1 + tr4;
	    ci2 = ti1 + ti4;
	    ci4 = ti1 - ti4;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 2] * cr2 - 
		    wa1[i - 1] * ci2;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 2] * ci2 + wa1[i 
		    - 1] * cr2;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 2] * cr3 - wa2[
		    i - 1] * ci3;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 2] * ci3 + wa2[i - 
		    1] * cr3;
	    ch[i - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 2] * cr4 - 
		    wa3[i - 1] * ci4;
	    ch[i + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 2] * ci4 + wa3[i 
		    - 1] * cr4;
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti1 = cc[((k << 2) + 2) * cc_dim1 + 1] + cc[((k << 2) + 4) * cc_dim1 
		+ 1];
	ti2 = cc[((k << 2) + 4) * cc_dim1 + 1] - cc[((k << 2) + 2) * cc_dim1 
		+ 1];
	tr1 = cc[*ido + ((k << 2) + 1) * cc_dim1] - cc[*ido + ((k << 2) + 3) *
		 cc_dim1];
	tr2 = cc[*ido + ((k << 2) + 1) * cc_dim1] + cc[*ido + ((k << 2) + 3) *
		 cc_dim1];
	ch[*ido + (k + ch_dim2) * ch_dim1] = tr2 + tr2;
	ch[*ido + (k + (ch_dim2 << 1)) * ch_dim1] = sqrt2 * (tr1 - ti1);
	ch[*ido + (k + ch_dim2 * 3) * ch_dim1] = ti2 + ti2;
	ch[*ido + (k + (ch_dim2 << 2)) * ch_dim1] = -(doublereal)sqrt2 * (tr1 
		+ ti1);
/* L106: */
    }
L107:
    return 0;
} /* radb4_ */

/* Subroutine */ int radb5_(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
integer *ido, *l1;
real *cc, *ch, *wa1, *wa2, *wa3, *wa4;
{
    /* Initialized data */

    static real tr11 = (float).309016994374947;
    static real ti11 = (float).951056516295154;
    static real tr12 = (float)-.809016994374947;
    static real ti12 = (float).587785252292473;

    /* System generated locals */
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i, k, ic;
    static real ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, 
	    ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
    static integer idp2;

    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 6 + 1;
    cc -= cc_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti5 = cc[(k * 5 + 3) * cc_dim1 + 1] + cc[(k * 5 + 3) * cc_dim1 + 1];
	ti4 = cc[(k * 5 + 5) * cc_dim1 + 1] + cc[(k * 5 + 5) * cc_dim1 + 1];
	tr2 = cc[*ido + (k * 5 + 2) * cc_dim1] + cc[*ido + (k * 5 + 2) * 
		cc_dim1];
	tr3 = cc[*ido + (k * 5 + 4) * cc_dim1] + cc[*ido + (k * 5 + 4) * 
		cc_dim1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 5 + 1) * cc_dim1 + 1] + tr2 
		+ tr3;
	cr2 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr11 * tr2 + tr12 * tr3;
	cr3 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr12 * tr2 + tr11 * tr3;
	ci5 = ti11 * ti5 + ti12 * ti4;
	ci4 = ti12 * ti5 - ti11 * ti4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci5;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr3 - ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = cr3 + ci4;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 1] = cr2 + ci5;
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    ti5 = cc[i + (k * 5 + 3) * cc_dim1] + cc[ic + (k * 5 + 2) * 
		    cc_dim1];
	    ti2 = cc[i + (k * 5 + 3) * cc_dim1] - cc[ic + (k * 5 + 2) * 
		    cc_dim1];
	    ti4 = cc[i + (k * 5 + 5) * cc_dim1] + cc[ic + (k * 5 + 4) * 
		    cc_dim1];
	    ti3 = cc[i + (k * 5 + 5) * cc_dim1] - cc[ic + (k * 5 + 4) * 
		    cc_dim1];
	    tr5 = cc[i - 1 + (k * 5 + 3) * cc_dim1] - cc[ic - 1 + (k * 5 + 2) 
		    * cc_dim1];
	    tr2 = cc[i - 1 + (k * 5 + 3) * cc_dim1] + cc[ic - 1 + (k * 5 + 2) 
		    * cc_dim1];
	    tr4 = cc[i - 1 + (k * 5 + 5) * cc_dim1] - cc[ic - 1 + (k * 5 + 4) 
		    * cc_dim1];
	    tr3 = cc[i - 1 + (k * 5 + 5) * cc_dim1] + cc[ic - 1 + (k * 5 + 4) 
		    * cc_dim1];
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = cc[i - 1 + (k * 5 + 1) * 
		    cc_dim1] + tr2 + tr3;
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * 5 + 1) * cc_dim1] + 
		    ti2 + ti3;
	    cr2 = cc[i - 1 + (k * 5 + 1) * cc_dim1] + tr11 * tr2 + tr12 * tr3;
	    ci2 = cc[i + (k * 5 + 1) * cc_dim1] + tr11 * ti2 + tr12 * ti3;
	    cr3 = cc[i - 1 + (k * 5 + 1) * cc_dim1] + tr12 * tr2 + tr11 * tr3;
	    ci3 = cc[i + (k * 5 + 1) * cc_dim1] + tr12 * ti2 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    dr3 = cr3 - ci4;
	    dr4 = cr3 + ci4;
	    di3 = ci3 + cr4;
	    di4 = ci3 - cr4;
	    dr5 = cr2 + ci5;
	    dr2 = cr2 - ci5;
	    di5 = ci2 - cr5;
	    di2 = ci2 + cr5;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 2] * dr2 - 
		    wa1[i - 1] * di2;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 2] * di2 + wa1[i 
		    - 1] * dr2;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 2] * dr3 - wa2[
		    i - 1] * di3;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 2] * di3 + wa2[i - 
		    1] * dr3;
	    ch[i - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 2] * dr4 - 
		    wa3[i - 1] * di4;
	    ch[i + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 2] * di4 + wa3[i 
		    - 1] * dr4;
	    ch[i - 1 + (k + ch_dim2 * 5) * ch_dim1] = wa4[i - 2] * dr5 - wa4[
		    i - 1] * di5;
	    ch[i + (k + ch_dim2 * 5) * ch_dim1] = wa4[i - 2] * di5 + wa4[i - 
		    1] * dr5;
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* radb5_ */

/* Subroutine */ int radbg_(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
integer *ido, *ip, *l1, *idl1;
real *cc, *c1, *c2, *ch, *ch2, *wa;
{
    /* Initialized data */

    static real tpi = (float)6.28318530717959;

    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_offset, c1_dim1,
	     c1_dim2, c1_offset, c2_dim1, c2_offset, ch2_dim1, ch2_offset, 
	    i__1, i__2, i__3;

    /* Builtin functions */
    double cos(), sin();

    /* Local variables */
    static integer idij, ipph, i, j, k, l, j2, ic, jc, lc, ik, is;
    static real dc2, ai1, ai2, ar1, ar2, ds2;
    static integer nbd;
    static real dcp, arg, dsp, ar1h, ar2h;
    static integer idp2, ipp2;

    /* Parameter adjustments */
    --wa;
    ch2_dim1 = *idl1;
    ch2_offset = ch2_dim1 + 1;
    ch2 -= ch2_offset;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    c2_dim1 = *idl1;
    c2_offset = c2_dim1 + 1;
    c2 -= c2_offset;
    c1_dim1 = *ido;
    c1_dim2 = *l1;
    c1_offset = c1_dim1 * (c1_dim2 + 1) + 1;
    c1 -= c1_offset;
    cc_dim1 = *ido;
    cc_dim2 = *ip;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;

    /* Function Body */
    arg = tpi / (real) (*ip);
    dcp = cos(arg);
    dsp = sin(arg);
    idp2 = *ido + 2;
    nbd = (*ido - 1) / 2;
    ipp2 = *ip + 2;
    ipph = (*ip + 1) / 2;
    if (*ido < *l1) {
	goto L103;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 1; i <= i__2; ++i) {
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L101: */
	}
/* L102: */
    }
    goto L106;
L103:
    i__1 = *ido;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L104: */
	}
/* L105: */
    }
L106:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = cc[*ido + (j2 - 2 + k * 
		    cc_dim2) * cc_dim1] + cc[*ido + (j2 - 2 + k * cc_dim2) * 
		    cc_dim1];
	    ch[(k + jc * ch_dim2) * ch_dim1 + 1] = cc[(j2 - 1 + k * cc_dim2) *
		     cc_dim1 + 1] + cc[(j2 - 1 + k * cc_dim2) * cc_dim1 + 1];
/* L107: */
	}
/* L108: */
    }
    if (*ido == 1) {
	goto L116;
    }
    if (nbd < *l1) {
	goto L112;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i = 3; i <= i__3; i += 2) {
		ic = idp2 - i;
		ch[i - 1 + (k + j * ch_dim2) * ch_dim1] = cc[i - 1 + ((j << 1)
			 - 1 + k * cc_dim2) * cc_dim1] + cc[ic - 1 + ((j << 1)
			 - 2 + k * cc_dim2) * cc_dim1];
		ch[i - 1 + (k + jc * ch_dim2) * ch_dim1] = cc[i - 1 + ((j << 
			1) - 1 + k * cc_dim2) * cc_dim1] - cc[ic - 1 + ((j << 
			1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i + (k + j * ch_dim2) * ch_dim1] = cc[i + ((j << 1) - 1 + 
			k * cc_dim2) * cc_dim1] - cc[ic + ((j << 1) - 2 + k * 
			cc_dim2) * cc_dim1];
		ch[i + (k + jc * ch_dim2) * ch_dim1] = cc[i + ((j << 1) - 1 + 
			k * cc_dim2) * cc_dim1] + cc[ic + ((j << 1) - 2 + k * 
			cc_dim2) * cc_dim1];
/* L109: */
	    }
/* L110: */
	}
/* L111: */
    }
    goto L116;
L112:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i - 1 + (k + j * ch_dim2) * ch_dim1] = cc[i - 1 + ((j << 1)
			 - 1 + k * cc_dim2) * cc_dim1] + cc[ic - 1 + ((j << 1)
			 - 2 + k * cc_dim2) * cc_dim1];
		ch[i - 1 + (k + jc * ch_dim2) * ch_dim1] = cc[i - 1 + ((j << 
			1) - 1 + k * cc_dim2) * cc_dim1] - cc[ic - 1 + ((j << 
			1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i + (k + j * ch_dim2) * ch_dim1] = cc[i + ((j << 1) - 1 + 
			k * cc_dim2) * cc_dim1] - cc[ic + ((j << 1) - 2 + k * 
			cc_dim2) * cc_dim1];
		ch[i + (k + jc * ch_dim2) * ch_dim1] = cc[i + ((j << 1) - 1 + 
			k * cc_dim2) * cc_dim1] + cc[ic + ((j << 1) - 2 + k * 
			cc_dim2) * cc_dim1];
/* L113: */
	    }
/* L114: */
	}
/* L115: */
    }
L116:
    ar1 = (float)1.;
    ai1 = (float)0.;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    c2[ik + l * c2_dim1] = ch2[ik + ch2_dim1] + ar1 * ch2[ik + (
		    ch2_dim1 << 1)];
	    c2[ik + lc * c2_dim1] = ai1 * ch2[ik + *ip * ch2_dim1];
/* L117: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		c2[ik + l * c2_dim1] += ar2 * ch2[ik + j * ch2_dim1];
		c2[ik + lc * c2_dim1] += ai2 * ch2[ik + jc * ch2_dim1];
/* L118: */
	    }
/* L119: */
	}
/* L120: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += ch2[ik + j * ch2_dim1];
/* L121: */
	}
/* L122: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * 
		    c1_dim1 + 1] - c1[(k + jc * c1_dim2) * c1_dim1 + 1];
	    ch[(k + jc * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * 
		    c1_dim1 + 1] + c1[(k + jc * c1_dim2) * c1_dim1 + 1];
/* L123: */
	}
/* L124: */
    }
    if (*ido == 1) {
	goto L132;
    }
    if (nbd < *l1) {
	goto L128;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i = 3; i <= i__3; i += 2) {
		ch[i - 1 + (k + j * ch_dim2) * ch_dim1] = c1[i - 1 + (k + j * 
			c1_dim2) * c1_dim1] - c1[i + (k + jc * c1_dim2) * 
			c1_dim1];
		ch[i - 1 + (k + jc * ch_dim2) * ch_dim1] = c1[i - 1 + (k + j *
			 c1_dim2) * c1_dim1] + c1[i + (k + jc * c1_dim2) * 
			c1_dim1];
		ch[i + (k + j * ch_dim2) * ch_dim1] = c1[i + (k + j * c1_dim2)
			 * c1_dim1] + c1[i - 1 + (k + jc * c1_dim2) * c1_dim1]
			;
		ch[i + (k + jc * ch_dim2) * ch_dim1] = c1[i + (k + j * 
			c1_dim2) * c1_dim1] - c1[i - 1 + (k + jc * c1_dim2) * 
			c1_dim1];
/* L125: */
	    }
/* L126: */
	}
/* L127: */
    }
    goto L132;
L128:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i - 1 + (k + j * ch_dim2) * ch_dim1] = c1[i - 1 + (k + j * 
			c1_dim2) * c1_dim1] - c1[i + (k + jc * c1_dim2) * 
			c1_dim1];
		ch[i - 1 + (k + jc * ch_dim2) * ch_dim1] = c1[i - 1 + (k + j *
			 c1_dim2) * c1_dim1] + c1[i + (k + jc * c1_dim2) * 
			c1_dim1];
		ch[i + (k + j * ch_dim2) * ch_dim1] = c1[i + (k + j * c1_dim2)
			 * c1_dim1] + c1[i - 1 + (k + jc * c1_dim2) * c1_dim1]
			;
		ch[i + (k + jc * ch_dim2) * ch_dim1] = c1[i + (k + j * 
			c1_dim2) * c1_dim1] - c1[i - 1 + (k + jc * c1_dim2) * 
			c1_dim1];
/* L129: */
	    }
/* L130: */
	}
/* L131: */
    }
L132:
    if (*ido == 1) {
	return 0;
    }
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L133: */
    }
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * 
		    ch_dim1 + 1];
/* L134: */
	}
/* L135: */
    }
    if (nbd > *l1) {
	goto L139;
    }
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i 
			- 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * ch[i 
			+ (k + j * ch_dim2) * ch_dim1];
		c1[i + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i + (
			k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i - 1 + (
			k + j * ch_dim2) * ch_dim1];
/* L136: */
	    }
/* L137: */
	}
/* L138: */
    }
    goto L143;
L139:
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i = 3; i <= i__3; i += 2) {
		idij += 2;
		c1[i - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i 
			- 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * ch[i 
			+ (k + j * ch_dim2) * ch_dim1];
		c1[i + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i + (
			k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i - 1 + (
			k + j * ch_dim2) * ch_dim1];
/* L140: */
	    }
/* L141: */
	}
/* L142: */
    }
L143:
    return 0;
} /* radbg_ */

/* Subroutine */ int radf2_(ido, l1, cc, ch, wa1)
integer *ido, *l1;
real *cc, *ch, *wa1;
{
    /* System generated locals */
    integer ch_dim1, ch_offset, cc_dim1, cc_dim2, cc_offset, i__1, i__2;

    /* Local variables */
    static integer i, k, ic;
    static real ti2, tr2;
    static integer idp2;

    /* Parameter adjustments */
    --wa1;
    ch_dim1 = *ido;
    ch_offset = ch_dim1 * 3 + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_dim2 = *l1;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[((k << 1) + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + 
		cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
	ch[*ido + ((k << 1) + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] 
		- cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    tr2 = wa1[i - 2] * cc[i - 1 + (k + (cc_dim2 << 1)) * cc_dim1] + 
		    wa1[i - 1] * cc[i + (k + (cc_dim2 << 1)) * cc_dim1];
	    ti2 = wa1[i - 2] * cc[i + (k + (cc_dim2 << 1)) * cc_dim1] - wa1[i 
		    - 1] * cc[i - 1 + (k + (cc_dim2 << 1)) * cc_dim1];
	    ch[i + ((k << 1) + 1) * ch_dim1] = cc[i + (k + cc_dim2) * cc_dim1]
		     + ti2;
	    ch[ic + ((k << 1) + 2) * ch_dim1] = ti2 - cc[i + (k + cc_dim2) * 
		    cc_dim1];
	    ch[i - 1 + ((k << 1) + 1) * ch_dim1] = cc[i - 1 + (k + cc_dim2) * 
		    cc_dim1] + tr2;
	    ch[ic - 1 + ((k << 1) + 2) * ch_dim1] = cc[i - 1 + (k + cc_dim2) *
		     cc_dim1] - tr2;
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[((k << 1) + 2) * ch_dim1 + 1] = -(doublereal)cc[*ido + (k + (
		cc_dim2 << 1)) * cc_dim1];
	ch[*ido + ((k << 1) + 1) * ch_dim1] = cc[*ido + (k + cc_dim2) * 
		cc_dim1];
/* L106: */
    }
L107:
    return 0;
} /* radf2_ */

/* Subroutine */ int radf3_(ido, l1, cc, ch, wa1, wa2)
integer *ido, *l1;
real *cc, *ch, *wa1, *wa2;
{
    /* Initialized data */

    static real taur = (float)-.5;
    static real taui = (float).866025403784439;

    /* System generated locals */
    integer ch_dim1, ch_offset, cc_dim1, cc_dim2, cc_offset, i__1, i__2;

    /* Local variables */
    static integer i, k, ic;
    static real ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;
    static integer idp2;

    /* Parameter adjustments */
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_offset = (ch_dim1 << 2) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_dim2 = *l1;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	cr2 = cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1];
	ch[(k * 3 + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + cr2;
	ch[(k * 3 + 3) * ch_dim1 + 1] = taui * (cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1] - cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1]);
	ch[*ido + (k * 3 + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + 
		taur * cr2;
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    dr2 = wa1[i - 2] * cc[i - 1 + (k + (cc_dim2 << 1)) * cc_dim1] + 
		    wa1[i - 1] * cc[i + (k + (cc_dim2 << 1)) * cc_dim1];
	    di2 = wa1[i - 2] * cc[i + (k + (cc_dim2 << 1)) * cc_dim1] - wa1[i 
		    - 1] * cc[i - 1 + (k + (cc_dim2 << 1)) * cc_dim1];
	    dr3 = wa2[i - 2] * cc[i - 1 + (k + cc_dim2 * 3) * cc_dim1] + wa2[
		    i - 1] * cc[i + (k + cc_dim2 * 3) * cc_dim1];
	    di3 = wa2[i - 2] * cc[i + (k + cc_dim2 * 3) * cc_dim1] - wa2[i - 
		    1] * cc[i - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    cr2 = dr2 + dr3;
	    ci2 = di2 + di3;
	    ch[i - 1 + (k * 3 + 1) * ch_dim1] = cc[i - 1 + (k + cc_dim2) * 
		    cc_dim1] + cr2;
	    ch[i + (k * 3 + 1) * ch_dim1] = cc[i + (k + cc_dim2) * cc_dim1] + 
		    ci2;
	    tr2 = cc[i - 1 + (k + cc_dim2) * cc_dim1] + taur * cr2;
	    ti2 = cc[i + (k + cc_dim2) * cc_dim1] + taur * ci2;
	    tr3 = taui * (di2 - di3);
	    ti3 = taui * (dr3 - dr2);
	    ch[i - 1 + (k * 3 + 3) * ch_dim1] = tr2 + tr3;
	    ch[ic - 1 + (k * 3 + 2) * ch_dim1] = tr2 - tr3;
	    ch[i + (k * 3 + 3) * ch_dim1] = ti2 + ti3;
	    ch[ic + (k * 3 + 2) * ch_dim1] = ti3 - ti2;
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* radf3_ */

/* Subroutine */ int radf4_(ido, l1, cc, ch, wa1, wa2, wa3)
integer *ido, *l1;
real *cc, *ch, *wa1, *wa2, *wa3;
{
    /* Initialized data */

    static real hsqt2 = (float).7071067811865475;

    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i, k, ic;
    static real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, 
	    tr3, tr4;
    static integer idp2;

    /* Parameter adjustments */
    --wa3;
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_offset = ch_dim1 * 5 + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_dim2 = *l1;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	tr1 = cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1] + cc[(k + (cc_dim2 << 2))
		 * cc_dim1 + 1];
	tr2 = cc[(k + cc_dim2) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1];
	ch[((k << 2) + 1) * ch_dim1 + 1] = tr1 + tr2;
	ch[*ido + ((k << 2) + 4) * ch_dim1] = tr2 - tr1;
	ch[*ido + ((k << 2) + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] 
		- cc[(k + cc_dim2 * 3) * cc_dim1 + 1];
	ch[((k << 2) + 3) * ch_dim1 + 1] = cc[(k + (cc_dim2 << 2)) * cc_dim1 
		+ 1] - cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    cr2 = wa1[i - 2] * cc[i - 1 + (k + (cc_dim2 << 1)) * cc_dim1] + 
		    wa1[i - 1] * cc[i + (k + (cc_dim2 << 1)) * cc_dim1];
	    ci2 = wa1[i - 2] * cc[i + (k + (cc_dim2 << 1)) * cc_dim1] - wa1[i 
		    - 1] * cc[i - 1 + (k + (cc_dim2 << 1)) * cc_dim1];
	    cr3 = wa2[i - 2] * cc[i - 1 + (k + cc_dim2 * 3) * cc_dim1] + wa2[
		    i - 1] * cc[i + (k + cc_dim2 * 3) * cc_dim1];
	    ci3 = wa2[i - 2] * cc[i + (k + cc_dim2 * 3) * cc_dim1] - wa2[i - 
		    1] * cc[i - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    cr4 = wa3[i - 2] * cc[i - 1 + (k + (cc_dim2 << 2)) * cc_dim1] + 
		    wa3[i - 1] * cc[i + (k + (cc_dim2 << 2)) * cc_dim1];
	    ci4 = wa3[i - 2] * cc[i + (k + (cc_dim2 << 2)) * cc_dim1] - wa3[i 
		    - 1] * cc[i - 1 + (k + (cc_dim2 << 2)) * cc_dim1];
	    tr1 = cr2 + cr4;
	    tr4 = cr4 - cr2;
	    ti1 = ci2 + ci4;
	    ti4 = ci2 - ci4;
	    ti2 = cc[i + (k + cc_dim2) * cc_dim1] + ci3;
	    ti3 = cc[i + (k + cc_dim2) * cc_dim1] - ci3;
	    tr2 = cc[i - 1 + (k + cc_dim2) * cc_dim1] + cr3;
	    tr3 = cc[i - 1 + (k + cc_dim2) * cc_dim1] - cr3;
	    ch[i - 1 + ((k << 2) + 1) * ch_dim1] = tr1 + tr2;
	    ch[ic - 1 + ((k << 2) + 4) * ch_dim1] = tr2 - tr1;
	    ch[i + ((k << 2) + 1) * ch_dim1] = ti1 + ti2;
	    ch[ic + ((k << 2) + 4) * ch_dim1] = ti1 - ti2;
	    ch[i - 1 + ((k << 2) + 3) * ch_dim1] = ti4 + tr3;
	    ch[ic - 1 + ((k << 2) + 2) * ch_dim1] = tr3 - ti4;
	    ch[i + ((k << 2) + 3) * ch_dim1] = tr4 + ti3;
	    ch[ic + ((k << 2) + 2) * ch_dim1] = tr4 - ti3;
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti1 = -(doublereal)hsqt2 * (cc[*ido + (k + (cc_dim2 << 1)) * cc_dim1] 
		+ cc[*ido + (k + (cc_dim2 << 2)) * cc_dim1]);
	tr1 = hsqt2 * (cc[*ido + (k + (cc_dim2 << 1)) * cc_dim1] - cc[*ido + (
		k + (cc_dim2 << 2)) * cc_dim1]);
	ch[*ido + ((k << 2) + 1) * ch_dim1] = tr1 + cc[*ido + (k + cc_dim2) * 
		cc_dim1];
	ch[*ido + ((k << 2) + 3) * ch_dim1] = cc[*ido + (k + cc_dim2) * 
		cc_dim1] - tr1;
	ch[((k << 2) + 2) * ch_dim1 + 1] = ti1 - cc[*ido + (k + cc_dim2 * 3) *
		 cc_dim1];
	ch[((k << 2) + 4) * ch_dim1 + 1] = ti1 + cc[*ido + (k + cc_dim2 * 3) *
		 cc_dim1];
/* L106: */
    }
L107:
    return 0;
} /* radf4_ */

/* Subroutine */ int radf5_(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
integer *ido, *l1;
real *cc, *ch, *wa1, *wa2, *wa3, *wa4;
{
    /* Initialized data */

    static real tr11 = (float).309016994374947;
    static real ti11 = (float).951056516295154;
    static real tr12 = (float)-.809016994374947;
    static real ti12 = (float).587785252292473;

    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i, k, ic;
    static real ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3, 
	    dr4, dr5, cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
    static integer idp2;

    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_offset = ch_dim1 * 6 + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_dim2 = *l1;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	cr2 = cc[(k + cc_dim2 * 5) * cc_dim1 + 1] + cc[(k + (cc_dim2 << 1)) * 
		cc_dim1 + 1];
	ci5 = cc[(k + cc_dim2 * 5) * cc_dim1 + 1] - cc[(k + (cc_dim2 << 1)) * 
		cc_dim1 + 1];
	cr3 = cc[(k + (cc_dim2 << 2)) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1];
	ci4 = cc[(k + (cc_dim2 << 2)) * cc_dim1 + 1] - cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1];
	ch[(k * 5 + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + cr2 
		+ cr3;
	ch[*ido + (k * 5 + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + 
		tr11 * cr2 + tr12 * cr3;
	ch[(k * 5 + 3) * ch_dim1 + 1] = ti11 * ci5 + ti12 * ci4;
	ch[*ido + (k * 5 + 4) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + 
		tr12 * cr2 + tr11 * cr3;
	ch[(k * 5 + 5) * ch_dim1 + 1] = ti12 * ci5 - ti11 * ci4;
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    dr2 = wa1[i - 2] * cc[i - 1 + (k + (cc_dim2 << 1)) * cc_dim1] + 
		    wa1[i - 1] * cc[i + (k + (cc_dim2 << 1)) * cc_dim1];
	    di2 = wa1[i - 2] * cc[i + (k + (cc_dim2 << 1)) * cc_dim1] - wa1[i 
		    - 1] * cc[i - 1 + (k + (cc_dim2 << 1)) * cc_dim1];
	    dr3 = wa2[i - 2] * cc[i - 1 + (k + cc_dim2 * 3) * cc_dim1] + wa2[
		    i - 1] * cc[i + (k + cc_dim2 * 3) * cc_dim1];
	    di3 = wa2[i - 2] * cc[i + (k + cc_dim2 * 3) * cc_dim1] - wa2[i - 
		    1] * cc[i - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    dr4 = wa3[i - 2] * cc[i - 1 + (k + (cc_dim2 << 2)) * cc_dim1] + 
		    wa3[i - 1] * cc[i + (k + (cc_dim2 << 2)) * cc_dim1];
	    di4 = wa3[i - 2] * cc[i + (k + (cc_dim2 << 2)) * cc_dim1] - wa3[i 
		    - 1] * cc[i - 1 + (k + (cc_dim2 << 2)) * cc_dim1];
	    dr5 = wa4[i - 2] * cc[i - 1 + (k + cc_dim2 * 5) * cc_dim1] + wa4[
		    i - 1] * cc[i + (k + cc_dim2 * 5) * cc_dim1];
	    di5 = wa4[i - 2] * cc[i + (k + cc_dim2 * 5) * cc_dim1] - wa4[i - 
		    1] * cc[i - 1 + (k + cc_dim2 * 5) * cc_dim1];
	    cr2 = dr2 + dr5;
	    ci5 = dr5 - dr2;
	    cr5 = di2 - di5;
	    ci2 = di2 + di5;
	    cr3 = dr3 + dr4;
	    ci4 = dr4 - dr3;
	    cr4 = di3 - di4;
	    ci3 = di3 + di4;
	    ch[i - 1 + (k * 5 + 1) * ch_dim1] = cc[i - 1 + (k + cc_dim2) * 
		    cc_dim1] + cr2 + cr3;
	    ch[i + (k * 5 + 1) * ch_dim1] = cc[i + (k + cc_dim2) * cc_dim1] + 
		    ci2 + ci3;
	    tr2 = cc[i - 1 + (k + cc_dim2) * cc_dim1] + tr11 * cr2 + tr12 * 
		    cr3;
	    ti2 = cc[i + (k + cc_dim2) * cc_dim1] + tr11 * ci2 + tr12 * ci3;
	    tr3 = cc[i - 1 + (k + cc_dim2) * cc_dim1] + tr12 * cr2 + tr11 * 
		    cr3;
	    ti3 = cc[i + (k + cc_dim2) * cc_dim1] + tr12 * ci2 + tr11 * ci3;
	    tr5 = ti11 * cr5 + ti12 * cr4;
	    ti5 = ti11 * ci5 + ti12 * ci4;
	    tr4 = ti12 * cr5 - ti11 * cr4;
	    ti4 = ti12 * ci5 - ti11 * ci4;
	    ch[i - 1 + (k * 5 + 3) * ch_dim1] = tr2 + tr5;
	    ch[ic - 1 + (k * 5 + 2) * ch_dim1] = tr2 - tr5;
	    ch[i + (k * 5 + 3) * ch_dim1] = ti2 + ti5;
	    ch[ic + (k * 5 + 2) * ch_dim1] = ti5 - ti2;
	    ch[i - 1 + (k * 5 + 5) * ch_dim1] = tr3 + tr4;
	    ch[ic - 1 + (k * 5 + 4) * ch_dim1] = tr3 - tr4;
	    ch[i + (k * 5 + 5) * ch_dim1] = ti3 + ti4;
	    ch[ic + (k * 5 + 4) * ch_dim1] = ti4 - ti3;
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* radf5_ */

/* Subroutine */ int radfg_(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
integer *ido, *ip, *l1, *idl1;
real *cc, *c1, *c2, *ch, *ch2, *wa;
{
    /* Initialized data */

    static real tpi = (float)6.28318530717959;

    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_offset, c1_dim1,
	     c1_dim2, c1_offset, c2_dim1, c2_offset, ch2_dim1, ch2_offset, 
	    i__1, i__2, i__3;

    /* Builtin functions */
    double cos(), sin();

    /* Local variables */
    static integer idij, ipph, i, j, k, l, j2, ic, jc, lc, ik, is;
    static real dc2, ai1, ai2, ar1, ar2, ds2;
    static integer nbd;
    static real dcp, arg, dsp, ar1h, ar2h;
    static integer idp2, ipp2;

    /* Parameter adjustments */
    --wa;
    ch2_dim1 = *idl1;
    ch2_offset = ch2_dim1 + 1;
    ch2 -= ch2_offset;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    c2_dim1 = *idl1;
    c2_offset = c2_dim1 + 1;
    c2 -= c2_offset;
    c1_dim1 = *ido;
    c1_dim2 = *l1;
    c1_offset = c1_dim1 * (c1_dim2 + 1) + 1;
    c1 -= c1_offset;
    cc_dim1 = *ido;
    cc_dim2 = *ip;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;

    /* Function Body */
    arg = tpi / (real) (*ip);
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (*ip + 1) / 2;
    ipp2 = *ip + 2;
    idp2 = *ido + 2;
    nbd = (*ido - 1) / 2;
    if (*ido == 1) {
	goto L119;
    }
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	ch2[ik + ch2_dim1] = c2[ik + c2_dim1];
/* L101: */
    }
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * 
		    c1_dim1 + 1];
/* L102: */
	}
/* L103: */
    }
    if (nbd > *l1) {
	goto L107;
    }
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i - 1 + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i 
			- 1 + (k + j * c1_dim2) * c1_dim1] + wa[idij] * c1[i 
			+ (k + j * c1_dim2) * c1_dim1];
		ch[i + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i + (
			k + j * c1_dim2) * c1_dim1] - wa[idij] * c1[i - 1 + (
			k + j * c1_dim2) * c1_dim1];
/* L104: */
	    }
/* L105: */
	}
/* L106: */
    }
    goto L111;
L107:
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i = 3; i <= i__3; i += 2) {
		idij += 2;
		ch[i - 1 + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i 
			- 1 + (k + j * c1_dim2) * c1_dim1] + wa[idij] * c1[i 
			+ (k + j * c1_dim2) * c1_dim1];
		ch[i + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i + (
			k + j * c1_dim2) * c1_dim1] - wa[idij] * c1[i - 1 + (
			k + j * c1_dim2) * c1_dim1];
/* L108: */
	    }
/* L109: */
	}
/* L110: */
    }
L111:
    if (nbd < *l1) {
	goto L115;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i = 3; i <= i__3; i += 2) {
		c1[i - 1 + (k + j * c1_dim2) * c1_dim1] = ch[i - 1 + (k + j * 
			ch_dim2) * ch_dim1] + ch[i - 1 + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i - 1 + (k + jc * c1_dim2) * c1_dim1] = ch[i + (k + j * 
			ch_dim2) * ch_dim1] - ch[i + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i + (k + j * c1_dim2) * c1_dim1] = ch[i + (k + j * ch_dim2)
			 * ch_dim1] + ch[i + (k + jc * ch_dim2) * ch_dim1];
		c1[i + (k + jc * c1_dim2) * c1_dim1] = ch[i - 1 + (k + jc * 
			ch_dim2) * ch_dim1] - ch[i - 1 + (k + j * ch_dim2) * 
			ch_dim1];
/* L112: */
	    }
/* L113: */
	}
/* L114: */
    }
    goto L121;
L115:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i - 1 + (k + j * c1_dim2) * c1_dim1] = ch[i - 1 + (k + j * 
			ch_dim2) * ch_dim1] + ch[i - 1 + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i - 1 + (k + jc * c1_dim2) * c1_dim1] = ch[i + (k + j * 
			ch_dim2) * ch_dim1] - ch[i + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i + (k + j * c1_dim2) * c1_dim1] = ch[i + (k + j * ch_dim2)
			 * ch_dim1] + ch[i + (k + jc * ch_dim2) * ch_dim1];
		c1[i + (k + jc * c1_dim2) * c1_dim1] = ch[i - 1 + (k + jc * 
			ch_dim2) * ch_dim1] - ch[i - 1 + (k + j * ch_dim2) * 
			ch_dim1];
/* L116: */
	    }
/* L117: */
	}
/* L118: */
    }
    goto L121;
L119:
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L120: */
    }
L121:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * 
		    ch_dim1 + 1] + ch[(k + jc * ch_dim2) * ch_dim1 + 1];
	    c1[(k + jc * c1_dim2) * c1_dim1 + 1] = ch[(k + jc * ch_dim2) * 
		    ch_dim1 + 1] - ch[(k + j * ch_dim2) * ch_dim1 + 1];
/* L122: */
	}
/* L123: */
    }

    ar1 = (float)1.;
    ai1 = (float)0.;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + l * ch2_dim1] = c2[ik + c2_dim1] + ar1 * c2[ik + (
		    c2_dim1 << 1)];
	    ch2[ik + lc * ch2_dim1] = ai1 * c2[ik + *ip * c2_dim1];
/* L124: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		ch2[ik + l * ch2_dim1] += ar2 * c2[ik + j * c2_dim1];
		ch2[ik + lc * ch2_dim1] += ai2 * c2[ik + jc * c2_dim1];
/* L125: */
	    }
/* L126: */
	}
/* L127: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += c2[ik + j * c2_dim1];
/* L128: */
	}
/* L129: */
    }

    if (*ido < *l1) {
	goto L132;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 1; i <= i__2; ++i) {
	    cc[i + (k * cc_dim2 + 1) * cc_dim1] = ch[i + (k + ch_dim2) * 
		    ch_dim1];
/* L130: */
	}
/* L131: */
    }
    goto L135;
L132:
    i__1 = *ido;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    cc[i + (k * cc_dim2 + 1) * cc_dim1] = ch[i + (k + ch_dim2) * 
		    ch_dim1];
/* L133: */
	}
/* L134: */
    }
L135:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    cc[*ido + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[(k + j * ch_dim2)
		     * ch_dim1 + 1];
	    cc[(j2 - 1 + k * cc_dim2) * cc_dim1 + 1] = ch[(k + jc * ch_dim2) *
		     ch_dim1 + 1];
/* L136: */
	}
/* L137: */
    }
    if (*ido == 1) {
	return 0;
    }
    if (nbd < *l1) {
	goto L141;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i = 3; i <= i__3; i += 2) {
		ic = idp2 - i;
		cc[i - 1 + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i - 1 + (k 
			+ j * ch_dim2) * ch_dim1] + ch[i - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[ic - 1 + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i - 1 + (k 
			+ j * ch_dim2) * ch_dim1] - ch[i - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[i + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i + (k + j * 
			ch_dim2) * ch_dim1] + ch[i + (k + jc * ch_dim2) * 
			ch_dim1];
		cc[ic + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i + (k + jc * 
			ch_dim2) * ch_dim1] - ch[i + (k + j * ch_dim2) * 
			ch_dim1];
/* L138: */
	    }
/* L139: */
	}
/* L140: */
    }
    return 0;
L141:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *ido;
	for (i = 3; i <= i__2; i += 2) {
	    ic = idp2 - i;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		cc[i - 1 + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i - 1 + (k 
			+ j * ch_dim2) * ch_dim1] + ch[i - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[ic - 1 + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i - 1 + (k 
			+ j * ch_dim2) * ch_dim1] - ch[i - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[i + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i + (k + j * 
			ch_dim2) * ch_dim1] + ch[i + (k + jc * ch_dim2) * 
			ch_dim1];
		cc[ic + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i + (k + jc * 
			ch_dim2) * ch_dim1] - ch[i + (k + j * ch_dim2) * 
			ch_dim1];
/* L142: */
	    }
/* L143: */
	}
/* L144: */
    }
    return 0;
} /* radfg_ */

