/********************************* ksync_rms.c *******************************
**
** Copyright 2001, University of Washington
**   This document contains private and confidential information and
**   its disclosure does not constitute publication.  All rights are
**   reserved by the University of Washington, the Baker Lab, 
**   Dylan Chivian, and Charlie Strauss, except those specifically 
**   granted by license.
**
**  Initial Author: Charlie Strauss (cems@lanl.gov)
**  $Revision: 1.1.1.1 $
**  $Date: 2001/11/28 03:33:02 $
**  $Author: dylan $
**
******************************************************************************/

/** Phil Bradley hacked the rms functions to return the rotation matrix, following **/
/**  along with 4rms.f in mammoth.   4/8/02 **/


/* rms.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <stdio.h>
#include <math.h>

#include "f2c.h"
/*  #include "ksync.h" */

/* Table of constant values */

#define PI 3.141592654

static double c_b5 = 1.;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;



double det3_(m)
double *m;
{
    /* System generated locals */
    double ret_val;

/* AUTHOR: charlie strauss 2001 */
/* determinant of a 3x3 matrix */
/* cems: cute factoid: det of a 3x3 is the dot product of one row with the cross product of the other two. */
/* this explains why a right hand coordinate system has a positive determinant. cute huh? */
    /* Parameter adjustments */
    m -= 4;

    /* Function Body */
    ret_val = m[10] * (m[5] * m[9] - m[8] * m[6]) - m[11] * (m[4] * m[9] - m[
	    7] * m[6]) + m[12] * (m[4] * m[8] - m[7] * m[5]);
    return ret_val;
} /* det3_ */



int rsym_evector(m,ev,mvec)
     double *m,*ev,mvec[3][3];
{
  static double xx,yy,xy,zx,yz,e1,e2,e3,znorm;
  int i;

  xy = m[3];
  zx = m[6];
  yz = m[7];
/*    xy = m[0][1]; */
/*    zx = m[0][2]; */
/*    yz = m[1][2]; */

  if ( ev[0] != ev[1] ) {
    for (i=0;i<2;i++) { 
      
      xx = m[0]-ev[i];
      yy = m[4]-ev[i];
/*        xx = m[0][0]-ev[i]; */
/*        yy = m[1][1]-ev[i]; */

      e1 = xy*yz-zx*yy;
      e2 = xy*zx-yz*xx;
      e3 = xx*yy-xy*xy;
            
      znorm= sqrt(e1*e1+e2*e2+e3*e3);
           
      mvec[0][i] = e1/znorm;
      mvec[1][i] = e2/znorm;
      mvec[2][i] = e3/znorm;
    }
    
    mvec[0][2] =  mvec[1][0]*mvec[2][1] -mvec[1][1]*mvec[2][0];
    mvec[1][2] = -mvec[0][0]*mvec[2][1] +mvec[0][1]*mvec[2][0];
    mvec[2][2] =  mvec[0][0]*mvec[1][1] -mvec[0][1]*mvec[1][0];
         
    return;
  }
  else {
    if (ev[1]!=ev[2]) {
      for (i=1;i<3;i++) {
	xx = m[0] - ev[i];
	yy = m[4] - ev[i];
/*  	xx = m[0][0] - ev[i]; */
/*  	yy = m[1][1] - ev[i]; */

	e1 = xy*yz-zx*yy;
	e2 = xy*zx-yz*xx;
	e3 = xx*yy-xy*xy;
               
	znorm= sqrt(e1*e1+e2*e2+e3*e3); 
               
	mvec[0][i] = e1/znorm;
	mvec[1][i] = e2/znorm;
	mvec[2][i] = e3/znorm;
      }
      mvec[0][0] =  mvec[1][1]*mvec[2][2] -mvec[1][2]*mvec[2][1];
      mvec[1][0] = -mvec[0][1]*mvec[2][2] +mvec[0][2]*mvec[2][1];
      mvec[2][0] =  mvec[0][1]*mvec[1][2] -mvec[0][2]*mvec[1][1];
      
      return;
    }
    else {
      for (i=0;i<3;i++) {
	mvec[0][i] = 0.0;
	mvec[1][i] = 0.0;
	mvec[2][i] = 0.0;
	mvec[i][i] = 0.0;
      }
      return;
    }
  }
}

int rsym_rotation(mm,m,ev,rot)
     
     double *mm,*m,*ev,rot[3][3];
{
  int i,j,k;
  static double temp[3][3],mvec[3][3],norm;

  rsym_evector(m,ev,mvec);

  for (i=0;i<2;i++) {
    norm = 1/sqrt(abs(ev[i]));

    for (j=0;j<3;j++) {
      temp[j][i] = 0.0;
      for (k=0;k<3;k++) temp[j][i] = temp[j][i]+mvec[k][i]*mm[j+3*k];
      temp[j][i] = temp[j][i]*norm;
    }
  }

  temp[0][2] =  temp[1][0]*temp[2][1] -temp[1][1]*temp[2][0];
  temp[1][2] = -temp[0][0]*temp[2][1] +temp[0][1]*temp[2][0];
  temp[2][2] =  temp[0][0]*temp[1][1] -temp[0][1]*temp[1][0];
  
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      rot[j][i] = 0.0;
      for (k=0;k<3;k++) {
	rot[j][i] = rot[j][i] + temp[i][k]*mvec[j][k];
      }
    }
  }
  return;
}




/* ======================================================================== */
/* Subroutine */ int rsym_eigenval__(m, ev)
double *m, *ev;
{
    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static double a, b, c;
    static double xx, xy, yy, xz, zz, yz;
    
    /* added for real cubic roots */
    double q, r, sq, theta;


/* computes the eigen values of a real symmetric 3x3 matrix */
/* AUTHOR: charlie strauss 2001 */
/* the method used is a deterministic analytic result that I hand factored.(whew!) */
/* Amusingly, while I'm suspect this factorization is not yet optimal in the number of calcs required */
/* I cannot find a signifcantly better one. */
/*  (if it were optimal I suspect I would not have to compute imaginary numbers that */
/*  I know must eventually cancel out to give a net real result.) */
/* this method relys on the fact that an analytic factoring of an order 3 polynomial exists. */
/* m(3,3) is a 3x3 real symmetric matrix: only the upper triangle is actually used */
/*  ev(3) is a real vector of eigen values, not neccesarily in sorted order. */
/* should be complex*16 but absoft sucks */
/* first, for sanity only, name some temporary variables */
    /* Parameter adjustments */
    --ev;
    m -= 4;

    /* Function Body */
    xx = m[4];
    yy = m[8];
    zz = m[12];
    xy = m[7];
    xz = m[10];
    yz = m[11];
/* coefficients of characterisitic polynomial */
    a = xx + yy + zz;
    b = -xx * zz - xx * yy - yy * zz + xy * xy + xz * xz + yz * yz;
    c = xx * yy * zz - xz * xz * yy - xy * xy * zz - yz * yz * xx + xy * 2 *
	     xz * yz;
/* eigenvals are the roots of the characteristic polymonial  0 = c+b*e +a*e^2  - e^3 */
/* solving for the three roots now: */
/* dont try to follow this in detail: its just a tricky */
/* factorization of the formulas for cubic equation roots. */

    /* real roots of cubic (sort of butt ugly) */
    b = -b;
    q = (a*a - 3.0*b) / 9.0;
    r = (2.0*a*a*a - 9.0*a*b + 27.0*c) / 54.0;
    if (q*q*q < r*r) {
/*        fprintf (stderr, "programmer f**ked up\n"); */
      return 1; /** catch this later on **/
    }
    sq = sqrt (q);
    theta = acos (r/(sq*sq*sq));
    sq *= 2.0;
    a /= -3.0;
    ev[1] = sq*cos(theta/3.0) -a;
    ev[2] = sq*cos((theta+2.0*PI)/3.0) -a;
    ev[3] = sq*cos((theta-2.0*PI)/3.0) -a;

    return 0;
} /* rsym_eigenval__ */



/********** RMS w/o SUPERPOSITION ********/
double rmsfit_(npoints, xx, yy)
integer *npoints;
double *xx, *yy;
{
    /* System generated locals */
    double ret_val;


    /* Local variables */
    static double det, xx_0__[15000]	/* was [3][5000] */, yy_0__[
	    15000]	/* was [3][5000] */;
    extern double rms2_();

    /* Fortran I/O blocks */
    static cilist io___44 = { 0, 0, 0, 0, 0 };


/*  this is just a wrapper for rms2() */
/*  its purpose is to hide the temporary arrays for COM offsets from the user. */
/*  unfortunately since were using fortran 77 and not fortran 90 we have to */
/*  allocate these temp arrays as big as we will ever need them */
/* real*8 ww(max_points) */
/* stuff that gets used by rms_setup */
/* outputs */
    /* Parameter adjustments */
    yy -= 4;
    xx -= 4;
    /* Function Body */
    if ( (*npoints) > 5000) {
      fprintf (stderr, "function: rmsfit: maximum number of points exceeded, %i; change max_points\n", (*npoints));
      exit (-2);
    }
    ret_val = rms2_(npoints, &xx[4], &yy[4], xx_0__, yy_0__, &det);
    return ret_val;
} /* rmsfit_ */





double rms2_(npoints, xx_0__, yy_0__, xx, yy, det)
integer *npoints;
double *xx_0__, *yy_0__, *xx, *yy, *det;
{
    /* System generated locals */
    integer i__1;
    double ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static double m_moment__[9]	/* was [3][3] */, temp1, temp2, temp3;
    static integer i__, j, k;
    static double rr_moment__[9]	/* was [3][3] */;
    extern /* Subroutine */ int rsym_eigenval__();
    static double ev[3], handedness, rms_ctx__;
    extern double det3_();
    static double rms_sum__;

/* AUTHOR: charlie strauss 2001 */
/*   computes the rms between two weighted point vectors. */
/*   xx_0,yy_0 are the input vectors of of points and ww is their weights */
/*   xx,yy are by-product output vectors of the same points offset to remove center of mass */
/*   det is an out value of the determinant of the cross moment matrix */
/*   returned value is the rms */

/*   most of this is double precision for good reasons.  first there are some large */
/*   differences of small numbers.  and second the rsymm_eignen() routine can internally have numbers */
/*   larger than the largest real*4 number.  (you could do some fancy foot work to rescale things if */
/*   you really had a problem with this. */

/*   NOTE: det is a double precision real */
/*   NOTE: (xx,yy) can be same arrays as (xx_0,yy_0) if desired */


/* ww(npoints) */
/* align center of mass to origin */
    /* Parameter adjustments */
    yy -= 4;
    xx -= 4;
    yy_0__ -= 4;
    xx_0__ -= 4;

    /* Function Body */
    for (k = 1; k <= 3; ++k) {
	temp1 = 0.;
	temp2 = 0.;
	temp3 = 0.;
	i__1 = *npoints;
	for (j = 1; j <= i__1; ++j) {
	    temp1 += xx_0__[k + j * 3];
/* *WW(j) */
	    temp2 += yy_0__[k + j * 3];
/* *WW(j) */
	    temp3 += 1.;
/* write(0,*) j,temp1,temp2,temp3 */
/* WW(j) */
/* L11: */
	}
	temp1 /= temp3;
	temp2 /= temp3;
	i__1 = *npoints;
	for (j = 1; j <= i__1; ++j) {
	    xx[k + j * 3] = xx_0__[k + j * 3] - temp1;
	    yy[k + j * 3] = yy_0__[k + j * 3] - temp2;
/* L10: */
	}
/* L12: */
    }
/* 	Make cross moments matrix   INCLUDE THE WEIGHTS HERE */
    for (k = 1; k <= 3; ++k) {
	for (j = 1; j <= 3; ++j) {
	    temp1 = (float)0.;
	    i__1 = *npoints;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		temp1 += yy[k + i__ * 3] * xx[j + i__ * 3];
/* *ww(i) */
	    }
	    m_moment__[k + j * 3 - 4] = temp1 / temp3;
/* rescale by temp3 (helps keeps numbers real*4 sized) */
	}
    }
/* write(0,*) ('m_m',m_moment(1,k),m_moment(2,k),m_moment(3,k),k=1,3) */
/* write(0,*) 'temp',temp1,temp2,temp3 */
    *det = det3_(m_moment__);
/* get handedness  of frame from determinant */
    if (abs(*det) <= (float)1e-24) {
/*  write(0,*) 'Warning:degenerate cross moments: det=',det */
/* might think about returning a zero rms, to avoid any chance of Floating Point Errors? */
	ret_val = 0.;
	return ret_val;
    }
    handedness = (*det < 0) ? -1.0 : 1.0;
/* weird but documented "feature" of DSIGN(a,b) (but not SIGN) is that if fails if a<0 */
/*  multiply cross moments by itself */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    rr_moment__[i__ + j * 3 - 4] = m_moment__[i__ * 3 - 3] * 
		    m_moment__[j * 3 - 3] + m_moment__[i__ * 3 - 2] * 
		    m_moment__[j * 3 - 2] + m_moment__[i__ * 3 - 1] * 
		    m_moment__[j * 3 - 1];
	}
    }
/*  compute eigen values of cross-cross moments */
    rsym_eigenval__(rr_moment__, ev);
    if (handedness < 0.) {
/* reorder eigen values  so that ev(1) is the smallest eigenvalue */
	if (ev[1] > ev[0]) {
	    if (ev[0] > ev[2]) {
		temp1 = ev[0];
		ev[0] = ev[2];
		ev[2] = temp1;
	    }
	} else {
	    if (ev[1] > ev[2]) {
		temp1 = ev[0];
		ev[0] = ev[2];
		ev[2] = temp1;
	    } else {
		temp1 = ev[0];
		ev[0] = ev[1];
		ev[1] = temp1;
	    }
	}
/* ev(1) is now the smallest eigen value.  the other two are not sorted. */
/* now we must catch the special case of the rotation with inversion. */
/* we cannot allow inversion rotations. */
/* fortunatley, and curiously, the optimal non-inverted rotation matrix */
/* will have the similar eigen values. */
/* we just have to make a slight change in how we handle things depending on determinant */
	
    }


    rms_ctx__ = sqrt((abs(ev[2]))) + sqrt((abs(ev[1]))) + handedness * sqrt((
	    abs(ev[0])));
/* the abs() are theoretically unneccessary since the eigen values of a real symmetric */
/* matrix are non-negative.  in practice sometimes small eigen vals end up just negative */
    rms_sum__ = 0.;
    i__1 = *npoints;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* Computing 2nd power */
	    d__1 = yy[j + i__ * 3];
/* Computing 2nd power */
	    d__2 = xx[j + i__ * 3];
	    rms_sum__ += d__1 * d__1 + d__2 * d__2;
/* *ww(i) */
	}
    }
    rms_sum__ /= temp3;
/* and combine the outer and cross terms into the final calculation. */
/*  (the abs() just saves us a headache when the roundoff error accidantally makes the sum negative) */
        
    ret_val = sqrt((d__1 = rms_sum__ - rms_ctx__ * 2., abs(d__1)));
    return ret_val;
} /* rms2_ */






/*********** SUPERPOSITIONS **************/

double fit_rmsfit_(npoints, xx, yy, R, T)

integer *npoints;
double *xx, *yy;
double R[3][3]; /** holds the rotation matrix **/
double T[3]; /** holds the translation: mapping of yy to xx-coord system is given by **/
           /** R*yy + T **/

{
    /* System generated locals */
    double ret_val;


    /* Local variables */
    static double det, xx_0__[15000]	/* was [3][5000] */, yy_0__[
	    15000]	/* was [3][5000] */;
    extern double fit_rms2_();

    /* Fortran I/O blocks */
    static cilist io___44 = { 0, 0, 0, 0, 0 };


/*  this is just a wrapper for rms2() */
/*  its purpose is to hide the temporary arrays for COM offsets from the user. */
/*  unfortunately since were using fortran 77 and not fortran 90 we have to */
/*  allocate these temp arrays as big as we will ever need them */
/* real*8 ww(max_points) */
/* stuff that gets used by rms_setup */
/* outputs */
    /* Parameter adjustments */
    yy -= 4;
    xx -= 4;
    /* Function Body */
    if ( (*npoints) > 5000) {
      fprintf (stderr, "function: rmsfit: maximum number of points exceeded; %i change max_points\n", *npoints);
      exit (-2);
    }
    ret_val = fit_rms2_(npoints, &xx[4], &yy[4], xx_0__, yy_0__, &det, R, T);
    return ret_val;
} /* rmsfit_ */




double fit_rms2_(npoints, xx_0__, yy_0__, xx, yy, det, R, T)
integer *npoints;
double *xx_0__, *yy_0__, *xx, *yy, *det;
double R[3][3]; /** holds the rotation matrix **/
double T[3]; /** holds the translation: mapping of yy to xx-coord system is given by **/
           /** R*yy + T **/
{
    /* System generated locals */
    integer i__1;
    double ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static double m_moment__[9]	/* was [3][3] */, temp1, temp2, temp3;
    static integer i__, j, k;
    static double rr_moment__[9]	/* was [3][3] */;
    extern /* Subroutine */ int rsym_eigenval__();
    static double ev[3], handedness, rms_ctx__;
    extern double det3_();
    static double rms_sum__;

    static double phil1,phil2,phil3, x_cm[3], y_cm[3];;

/* AUTHOR: charlie strauss 2001 */
/*   computes the rms between two weighted point vectors. */
/*   xx_0,yy_0 are the input vectors of of points and ww is their weights */
/*   xx,yy are by-product output vectors of the same points offset to remove center of mass */
/*   det is an out value of the determinant of the cross moment matrix */
/*   returned value is the rms */

/*   most of this is double precision for good reasons.  first there are some large */
/*   differences of small numbers.  and second the rsymm_eignen() routine can internally have numbers */
/*   larger than the largest real*4 number.  (you could do some fancy foot work to rescale things if */
/*   you really had a problem with this. */

/*   NOTE: det is a double precision real */
/*   NOTE: (xx,yy) can be same arrays as (xx_0,yy_0) if desired */


/* ww(npoints) */
/* align center of mass to origin */
    /* Parameter adjustments */
    yy -= 4;
    xx -= 4;
    yy_0__ -= 4;
    xx_0__ -= 4;

    /* Function Body */
    for (k = 1; k <= 3; ++k) {
	temp1 = 0.;
	temp2 = 0.;
	temp3 = 0.;
	i__1 = *npoints;
	for (j = 1; j <= i__1; ++j) {
	    temp1 += xx_0__[k + j * 3];
/* *WW(j) */
	    temp2 += yy_0__[k + j * 3];
/* *WW(j) */
	    temp3 += 1.;
/* write(0,*) j,temp1,temp2,temp3 */
/* WW(j) */
/* L11: */
	}
	temp1 /= temp3;
	temp2 /= temp3;
	i__1 = *npoints;
	for (j = 1; j <= i__1; ++j) {
	    xx[k + j * 3] = xx_0__[k + j * 3] - temp1;
	    yy[k + j * 3] = yy_0__[k + j * 3] - temp2;
/* L10: */
	}
/* L12: */
	x_cm[k-1] = temp1;
	y_cm[k-1] = temp2;
    }
/* 	Make cross moments matrix   INCLUDE THE WEIGHTS HERE */
    for (k = 1; k <= 3; ++k) {
	for (j = 1; j <= 3; ++j) {
	    temp1 = (float)0.;
	    i__1 = *npoints;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		temp1 += yy[k + i__ * 3] * xx[j + i__ * 3];
/* *ww(i) */
	    }
	    m_moment__[k + j * 3 - 4] = temp1 / temp3;
/* rescale by temp3 (helps keeps numbers real*4 sized) */
	}
    }
/* write(0,*) ('m_m',m_moment(1,k),m_moment(2,k),m_moment(3,k),k=1,3) */
/* write(0,*) 'temp',temp1,temp2,temp3 */
    *det = det3_(m_moment__);
/* get handedness  of frame from determinant */
    if (abs(*det) <= (float)1e-24) {
/*  write(0,*) 'Warning:degenerate cross moments: det=',det */
/* might think about returning a zero rms, to avoid any chance of Floating Point Errors? */
      ret_val = 2000.0;
      return ret_val;
    }
    handedness = (*det < 0) ? -1.0 : 1.0;
/* weird but documented "feature" of DSIGN(a,b) (but not SIGN) is that if fails if a<0 */
/*  multiply cross moments by itself */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    rr_moment__[i__ + j * 3 - 4] = m_moment__[i__ * 3 - 3] * 
		    m_moment__[j * 3 - 3] + m_moment__[i__ * 3 - 2] * 
		    m_moment__[j * 3 - 2] + m_moment__[i__ * 3 - 1] * 
		    m_moment__[j * 3 - 1];
	}
    }


/*  compute eigen values of cross-cross moments */
    
    if (  rsym_eigenval__(rr_moment__, ev)  ) { /** should return 0 **/
      return 1000.0;
    }
    

/*      if (handedness < 0.) { */
    if (1) {
      
/* reorder eigen values  so that ev(3!!!) is the smallest eigenvalue */
        if (ev[1] > ev[2]) {
	    if (ev[2] > ev[0]) {
		temp1 = ev[2];
		ev[2] = ev[0];
		ev[0] = temp1;
	    }
	} else {
	    if (ev[1] > ev[0]) {
		temp1 = ev[2];
		ev[2] = ev[0];
		ev[0] = temp1;
	    } else {
		temp1 = ev[2];
		ev[2] = ev[1];
		ev[1] = temp1;
	    }
	}
/* ev(3) is now the smallest eigen value.  the other two are not sorted. */
/* now we must catch the special case of the rotation with inversion. */
/* we cannot allow inversion rotations. */
/* fortunatley, and curiously, the optimal non-inverted rotation matrix */
/* will have the similar eigen values. */
/* we just have to make a slight change in how we handle things depending on determinant */
    }


    /** Here we do rsym_rotation(m_moment,rr_moment,ev,r) **/ 


    rsym_rotation(m_moment__,rr_moment__,ev,R);

    phil1 = R[0][0]*y_cm[0] + R[0][1]*y_cm[1] + R[0][2]*y_cm[2];
    phil2 = R[1][0]*y_cm[0] + R[1][1]*y_cm[1] + R[1][2]*y_cm[2];
    phil3 = R[2][0]*y_cm[0] + R[2][1]*y_cm[1] + R[2][2]*y_cm[2];
    T[0] = x_cm[0] - phil1;
    T[1] = x_cm[1] - phil2;
    T[2] = x_cm[2] - phil3;
    




    rms_ctx__ = sqrt((abs(ev[0]))) + sqrt((abs(ev[1]))) + handedness * sqrt((
	    abs(ev[2])));
/* the abs() are theoretically unneccessary since the eigen values of a real symmetric */
/* matrix are non-negative.  in practice sometimes small eigen vals end up just negative */
    rms_sum__ = 0.;
    i__1 = *npoints;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* Computing 2nd power */
	    d__1 = yy[j + i__ * 3];
/* Computing 2nd power */
	    d__2 = xx[j + i__ * 3];
	    rms_sum__ += d__1 * d__1 + d__2 * d__2;
/* *ww(i) */
	}
    }
    rms_sum__ /= temp3;
/* and combine the outer and cross terms into the final calculation. */
/*  (the abs() just saves us a headache when the roundoff error accidantally makes the sum negative) */
    ret_val = sqrt((d__1 = rms_sum__ - rms_ctx__ * 2., abs(d__1)));
    return ret_val;
} /* rms2_ */

