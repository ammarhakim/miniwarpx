/* tm_sys_eig.f -- translated by f2c (version 20060506).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__100 = 100;
static doublereal c_b7 = 1.;
static doublereal c_b8 = 0.;


/*     ================================================================== */
/* Subroutine */ int tm_sys_eig__(integer *meqn, doublereal *qli, doublereal *
	qri, doublereal *wr, doublereal *vl, doublereal *vr)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal q[10], r__, u, v, w, wi[10], vll[100]	/* was [10][
	    10] */, wxx, wxy, wxz, wyy, wyz, wzz, dfdq[100]	/* was [10][
	    10] */;
    static integer info;
    static doublereal work[1000], vlvr[100]	/* was [10][10] */;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgeev_(char *, char *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);

/*     ================================================================== */

/*     Computes flux eigensystem for 10 moment equations */

/*     INTENT(IN) */
/*     INTENT(OUT) */
/*     local variables/arrays */
/*     stuff needed for LAPACK routines */
/*      integer iwork,ipiv(meqn) */

/*     iwork = meqn*100          ! LAPACK work space size */
/*     compute average state at interface */
    /* Parameter adjustments */
    vr -= 11;
    vl -= 11;
    --wr;
    --qri;
    --qli;

    /* Function Body */
    i__1 = *meqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ - 1] = (qli[i__] + qri[i__]) * .5;
    }
/*     compute primitive variables */
    r__ = q[0];
    u = q[1] / q[0];
    v = q[2] / q[0];
    w = q[3] / q[0];
    wxx = (q[4] - r__ * u * u) / r__;
    wxy = (q[5] - r__ * u * v) / r__;
    wxz = (q[6] - r__ * u * w) / r__;
    wyy = (q[7] - r__ * v * v) / r__;
    wyz = (q[8] - r__ * v * w) / r__;
    wzz = (q[9] - r__ * w * w) / r__;
/*     compute flux jacobian */
    i__1 = *meqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *meqn;
	for (j = 1; j <= i__2; ++j) {
/*     ----- initialize everything to zero */
	    dfdq[i__ + j * 10 - 11] = 0.;
	}
    }
/*     ++ row 1: */
    dfdq[10] = 1.;
/*     ++ row 2 */
    dfdq[41] = 1.;
/*     ++ row 3 */
    dfdq[52] = 1.;
/*     ++ row 4 */
    dfdq[63] = 1.;
/*     ++ row 5 */
/* Computing 3rd power */
    d__1 = u;
    dfdq[4] = d__1 * (d__1 * d__1) - u * 3. * wxx;
/* Computing 2nd power */
    d__1 = u;
    dfdq[14] = d__1 * d__1 * -3. + wxx * 3.;
    dfdq[44] = u * 3.;
/*     ++ row 6 */
/* Computing 2nd power */
    d__1 = u;
    dfdq[5] = v * (d__1 * d__1) - u * 2. * wxy - v * wxx;
    dfdq[15] = u * -2. * v + wxy * 2.;
/* Computing 2nd power */
    d__1 = u;
    dfdq[25] = -(d__1 * d__1) + wxx;
    dfdq[45] = v;
    dfdq[55] = u * 2.;
/*     ++ row 7 */
/* Computing 2nd power */
    d__1 = u;
    dfdq[6] = w * (d__1 * d__1) - u * 2. * wxz - w * wxx;
    dfdq[16] = u * -2. * w + wxz * 2.;
/* Computing 2nd power */
    d__1 = u;
    dfdq[36] = -(d__1 * d__1) + wxx;
    dfdq[46] = w;
    dfdq[66] = u * 2.;
/*     ++ row 8 */
/* Computing 2nd power */
    d__1 = v;
    dfdq[7] = u * (d__1 * d__1) - u * wyy - v * 2. * wxy;
/* Computing 2nd power */
    d__1 = v;
    dfdq[17] = -(d__1 * d__1) + wyy;
    dfdq[27] = u * -2. * v + wxy * 2.;
    dfdq[57] = v * 2.;
    dfdq[77] = u;
/*     ++ row 9 */
    dfdq[8] = u * v * w - u * wyz - v * wxz - w * wxy;
    dfdq[18] = -v * w + wyz;
    dfdq[28] = -u * w + wxz;
    dfdq[38] = -u * v + wxy;
    dfdq[58] = w;
    dfdq[68] = v;
    dfdq[88] = u;
/*     ++ row 10 */
/* Computing 2nd power */
    d__1 = w;
    dfdq[9] = u * (d__1 * d__1) - u * wzz - w * 2. * wxz;
/* Computing 2nd power */
    d__1 = w;
    dfdq[19] = -(d__1 * d__1) + wzz;
    dfdq[39] = u * -2. * w + wxz * 2.;
    dfdq[69] = w * 2.;
    dfdq[99] = u;
/*     compute eigenvalues and right eigenvectors */
    dgeev_("V", "V", meqn, dfdq, meqn, &wr[1], wi, vll, meqn, &vr[11], meqn, 
	    work, &c__100, &info, (ftnlen)1, (ftnlen)1);
/*     find the product of the left and right eigenvectors. This should */
/*     be a diagonal matrix */
    dgemm_("T", "N", meqn, meqn, meqn, &c_b7, vll, meqn, &vr[11], meqn, &c_b8,
	     vlvr, meqn, (ftnlen)1, (ftnlen)1);
/*     renormalize left eigenvectors to give correct normalization */
/*     (vl'*vr = unit matrix) */
    i__1 = *meqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*         write(*,*) i, vlvr(i,i) */
	i__2 = *meqn;
	for (j = 1; j <= i__2; ++j) {
	    vll[i__ + j * 10 - 11] /= vlvr[j + j * 10 - 11];
	}
    }
/*     transpose left eigenvector matrix, */
    i__1 = *meqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *meqn;
	for (j = 1; j <= i__2; ++j) {
	    vl[i__ + j * 10] = vll[j + i__ * 10 - 11];
	}
    }
    return 0;
/*     -- subroutine tm_sys_eig ------------------------------------------ */
} /* tm_sys_eig__ */

