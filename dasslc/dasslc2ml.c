/*
 * $Log:        dasslc2ml.c,v $
 * Revision 1.0  07/10/16  20:05 {arge}
 * DASSLC version
 *
 * Revision 1.0  07/10/16  20:05 {arge}
 * - created mex gateway to dasslc DAE solver
 * Revision 1.1  08/04/29  10:15 {arge}
 * - setting up to compile for MATLAB versions > 5.3
 * Revision 1.2  08/05/18  11:00 {arge}
 * - adding optional argument yp0 (initial yp)
 */

/*
  **************************************************************************
  * INTERFACE FOR MATLAB AND DASSLC (DIFFERENTIAL-ALGEBRAIC SYSTEM SOLVER) *
  **************************************************************************

			  (c) Copyright 2007, 2008

			  Argimiro R. Secchi, UFRGS

			Simulation Laboratory - LASIM
            GIMSCOP (Group of Integration, Modeling, Simulation,
			         Control, and Optimization of Processes)
		      Department of Chemical Engineering
		    Federal University of Rio Grande do Sul
		      Porto Alegre, RS - Brazil
			 e-mail: arge@enq.ufrgs.br

  This code is published under the BSD License
  http://www.opensource.org/licenses/bsd-license.php
*/

#include <io.h>
#include <stdio.h>
#include <string.h>
#include "dasslc.h"
#undef REAL
#include "mex.h"

#ifdef V4_COMPAT
# define mxArray Matrix
# define mxReal REAL
#endif

FILE *fout, *ferr;
#define CLOSEALL() (_flushall(), fclose(ferr), fclose(fout))

DASSLC_RES residuals;
DASSLC_JAC jacobian;

#define Y(i,j)	(*(py+neq*(i)+j))
#define YP(i,j)	(*(pyp+neq*(i)+j))

typedef struct userdata
{
	char *fname;
	char *jac;
	mxArray *prhs[5], *pjac, *prpar;
	int mpar, npar;
	double *y, *yp, *rpar;
} UDATA;

/*
 mexFunction: entry point for the gateway to MATLAB
 
 usage:
 [t,y,yp,output] = dasslc(f,tspan,y0,yp0,rpar,rtol,atol,index,inputfile,jac);

*/

#ifdef __STDC__
void mexFunction(
	int		nlhs,
	mxArray	*plhs[],
	int		nrhs,
#if !defined(V4_COMPAT)
    const mxArray *prhs[]         /* array of pointers to input arguments */
#else
    mxArray *prhs[]         /* array of pointers to input arguments */
#endif
	)
#else
void mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
mxArray *plhs[];
#if !defined(V4_COMPAT)
    const mxArray *prhs[];         /* array of pointers to input arguments */
#else
    mxArray *prhs[];         /* array of pointers to input arguments */
#endif
#endif
{
 PTR_ROOT root;
 double	t, tout, tf, *tspan, *y0, *yp0=NULL, *rpar=NULL, rtol=-1, atol=-1, *y, *yp, *pt, *py, *pyp;
 FAST double *ptr;
 int nt, nout, neq, m, n, output, *index=NULL;
 FAST int i, j;
 char fname[128], inputfile[128]={'?','\0'}, jac[128]={'\0'};
 UDATA user;

 if (nrhs < 3) mexErrMsgTxt("dasslc: minimum number of inputs is 3: (f,tspan,y0)!");
 if (nrhs > 10) mexErrMsgTxt("dasslc: minimum number of inputs is 10: (f,tspan,y0,yp0,rpar,rtol,atol,index,inputfile,jac)!");
 if (nlhs < 2) mexErrMsgTxt("dasslc: minimum number of outputs is 2: [t,y]");
 if (nlhs > 4) mexErrMsgTxt("dasslc: maximum number of outputs is 4: [t,y,yp,output]");

 m = mxGetM(prhs[2]);
 n = mxGetN(prhs[2]);
 neq = MAX(m,n); // number of variables
 if (MIN(m,n) != 1)	mexErrMsgTxt("dasslc: initial condition must be a vector, not an array!");

 m = mxGetM(prhs[1]);
 n = mxGetN(prhs[1]);
 nt = MAX(m,n); // number of timepoints
 if (MIN(m,n) != 1)	mexErrMsgTxt("dasslc: timespan must be a vector, not an array!");

  /* Assign pointers to the various parameters */
 switch (nrhs)
	{
	 case 10: if (!mxIsEmpty(prhs[9]))
			   {
				mxGetString(prhs[9], jac, 127);
				user.pjac = mxCreateDoubleMatrix(1, 1, mxREAL); // cj
			   }
	 case 9: if (!mxIsEmpty(prhs[8])) mxGetString(prhs[8], inputfile, 127);
 	 case 8: if (!mxIsEmpty(prhs[7]))
			   {
				ptr = (double *)mxGetPr(prhs[7]);
				m = mxGetM(prhs[7]);
				n = mxGetN(prhs[7]);
				if (MAX(m,n) != neq || MIN(m,n) != 1)
				  mexErrMsgTxt("dasslc: index of variables must be of same size of y0!");

				index = NEW(int, neq, "index");
				for (i = 0; i < neq; i++) index[i] = (int)ptr[i];
			   }
	 case 7: if (!mxIsEmpty(prhs[6])) atol = (double)*mxGetPr(prhs[6]);
	 case 6: if (!mxIsEmpty(prhs[5])) rtol = (double)*mxGetPr(prhs[5]);
	 case 5: if (!mxIsEmpty(prhs[4]))
			   {
				rpar = (double *)mxGetPr(prhs[4]);
				user.mpar = mxGetM(prhs[4]);
				user.npar = mxGetN(prhs[4]);
			   }
	 case 4: if (!mxIsEmpty(prhs[3])) yp0 = (double *)mxGetPr(prhs[3]);
	 default:
			 y0 = (double *)mxGetPr(prhs[2]);
			 tspan = (double *)mxGetPr(prhs[1]);
			 mxGetString(prhs[0], fname, 127);
	}

 if (!(fout = fopen("dasslc_out.txt","w"), fout)) mexErrMsgTxt("dasslc: could not open output file!");
 if (!(ferr = fopen("dasslc_err.txt","w"), fout)) {fclose(fout); mexErrMsgTxt("dasslc: could not open error file!");}

 user.fname = fname;
 user.jac = jac;
 user.rpar = rpar;
 root.user = (void *)&user; // residuals and jacobian functions names

 if (nt == 1)
   {
    t = 0;
	tf = tspan[0];
	tout = tf;
   }
 else
   {
    t = tspan[0];
    tf = tspan[nt-1];
	tout = tspan[1];
   }

 user.prhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); // t
 user.prhs[1] = mxCreateDoubleMatrix(neq, 1, mxREAL); // y
 user.prhs[2] = mxCreateDoubleMatrix(neq, 1, mxREAL); // yp
 if (user.rpar)	user.prpar = prhs[4];

 yp = NEW(double, neq, "yp");
 y = NEW(double, neq, "y");

 user.y = (double *)mxGetPr(user.prhs[1]);
 user.yp = (double *)mxGetPr(user.prhs[2]);

 memcpy((void *)y, (void *)y0, neq * sizeof(double));

 if (yp0) memcpy((void *)yp, (void *)yp0, neq * sizeof(double));
 else for (j = 0; j < neq; j++) yp[j] = 0; // initial guess

 if (output = daSetup (inputfile, &root, residuals, neq, t, y, yp, index, NULL, NULL), output)
   {CLOSEALL(); mexErrMsgTxt("dasslc: error in Setup stage!");}
 
 if (nrhs > 5)
	{
	 if (rtol >= 0)
		{
		 root.iter.rtol[0] = rtol;
		 root.iter.stol = 1;
		}
	 if (nrhs > 6 && atol >= 0)
	   {
		root.iter.atol[0] = atol;
		root.iter.stol = 1;
	   }
	}

 if (!yp0 && (output = dasslc (INITIAL_COND, &root, residuals, &t, tout, (*jac ? jacobian : NULL), NULL), output < 0))
   {CLOSEALL(); mexErrMsgTxt("dasslc: error in finding consistent initial condition!");}

 root.iter.istall = 0;  // no intermediate time steps

 if (nt == 1 || (nt == 2 && tout <= 1))
   {
    nout = 100;
	root.iter.istall = 1;  // intermediate time steps
   }
 else if (nt == 2) nout = (int)(tout + 0.5)+1;
 else nout = nt;

 pt = NEW(double, nout, "pt");
 *pt = t;

 py = NEW(double, neq*nout, "py");
 memcpy((void *)py, (void *)y, neq * sizeof(double));

 if (nlhs > 2)
   {
	pyp = NEW(double, neq*nout, "pyp");
    memcpy((void *)pyp, (void *)yp, neq * sizeof(double));
   }

 for (i = 1; tout <= tf;)
    {
	 if (output = dasslc (TRANSIENT, &root, residuals, &t, tout, (*jac ? jacobian : NULL), NULL), output < 0)
	   {
		CLOSEALL(); mexErrMsgTxt("dasslc: error during integration!");
		break;
	   }

     pt[i] = t;
	 memcpy((void *)(py+neq*i), (void *)y, neq * sizeof(double));
	 if (nlhs > 2) memcpy((void *)(pyp+neq*i), (void *)yp, neq * sizeof(double));

	 i++;
	 if (t == tout)
	   {
		if (nt == nout && tout < tf) tout += tspan[i] - tspan[i-1];
		else tout += 0.1; // finished intermediate timesteps
	   }

     if (i == nout && tout <= tf)
	   {
	    nout += 100;
		pt = RENEW(double, pt, nout, "t");
		py = RENEW(double, py, neq*nout, "y");
		if (nlhs > 2) pyp = RENEW(double, pyp, neq*nout, "yp");
	   }
	}

 nout = i;

 switch (nlhs)
	{
	 case 4: plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
			 *(mxGetPr(plhs[3])) = output;
	 case 3: plhs[2] = mxCreateDoubleMatrix(nout, neq, mxREAL);
			 ptr = (double *)mxGetPr(plhs[2]);
			 for (j = 0; j < neq; j++) for (i = 0; i < nout; i++, ptr++) *ptr = YP(i,j); 

	 default: plhs[1] = mxCreateDoubleMatrix(nout, neq, mxREAL);
			  ptr = (double *)mxGetPr(plhs[1]);
			  for (j = 0; j < neq; j++) for (i = 0; i < nout; i++, ptr++) *ptr = Y(i,j);
			  plhs[0] = mxCreateDoubleMatrix(nout, 1, mxREAL);
		      memcpy((void *)mxGetPr(plhs[0]), (void *)pt, nout * sizeof(double));
	}

 free(y);
 free(yp);
 free(pt);
 free(py);
 if (nlhs > 2) free(pyp);
 if (index) free(index);

 mxDestroyArray(user.prhs[0]);
 mxDestroyArray(user.prhs[1]);
 mxDestroyArray(user.prhs[2]);
 if (*jac) mxDestroyArray(user.pjac);

 if (root.filename) daStat (root.savefile, &root);

 daFree (&root);

 CLOSEALL();
}

/*
 Initialization function
*/
BOOL init_dasslc (PTR_ROOT *root)
{
	root -> iter.stol = 1; // scalar tolerance
	root -> filename = NULL; // no savefile

	return STAT_OK;
}

/*
 Calling function to MATLAB residuals function (F(t,y,yp)=0).
 
 usage: function [res,ires]=fun(t,y,yp,rpar)
*/
BOOL
residuals (PTR_ROOT *root, double t, double *y, double *yp, double *res, BOOL *jac)
{
 BOOL error;
 int m, n;
 mxArray *plhs[2];
 UDATA *user = (UDATA *)root -> user;

 *(mxGetPr(user -> prhs[0])) = t;
 memcpy((void *)user -> y, (void *)y, root -> rank * sizeof(double));
 memcpy((void *)user -> yp, (void *)yp, root -> rank * sizeof(double));

 if (!user -> rpar) n = 3;
 else
   {
    n = 4;
	user -> prhs[3] = user -> prpar;
   }

 mexCallMATLAB(2, plhs, n, user -> prhs, user -> fname);

 error = (BOOL)*mxGetPr(plhs[1]);
 m = mxGetM(plhs[0]);
 n = mxGetN(plhs[0]);

 if (MAX(m,n) != root->rank || MIN(m,n) != 1)
   mexErrMsgTxt("dasslc: residuals require a vector of the same size of variables!");

 memcpy((void *)res, (void *)mxGetPr(plhs[0]), root -> rank * sizeof(double));

 mxDestroyArray(plhs[0]);
 mxDestroyArray(plhs[1]);

 return error;
}					/* residuals */

/*
 Calling function to MATLAB jacobian function (cj*dF_dyp(t,y,yp)+dF_dy). *** FULL MATRIX ONLY YET ***
 
 usage: function [jac,ires]=fun(t,y,yp,cj,rpar)
		where jac is the transpose of the iteration matrix
*/
BOOL
jacobian (PTR_ROOT *root, double t, double *y, double *yp, double cj, void *ja, DASSLC_RES *residuals)
{
 BOOL error;
 int m, n;
 mxArray *plhs[2];
 UDATA *user = (UDATA *)root -> user;

 *(mxGetPr(user -> prhs[0])) = t;
 memcpy((void *)user -> y, (void *)y, root -> rank * sizeof(double));
 memcpy((void *)user -> yp, (void *)yp, root -> rank * sizeof(double));
 user -> prhs[3] = user -> pjac;
 *(mxGetPr(user -> pjac)) = cj;

 if (!user -> rpar) n = 4;
 else
   {
    n = 5;
	user -> prhs[4] = user -> prpar;
   }

 mexCallMATLAB(2, plhs, n, user -> prhs, user -> jac);

 error = (BOOL)*mxGetPr(plhs[1]);
 m = mxGetM(plhs[0]);
 n = mxGetN(plhs[0]);

 if (m != root->rank || n != root->rank)
   mexErrMsgTxt("dasslc: jacobian require a square matrix on the size of variables!");

 memcpy((void *)ja, (void *)mxGetPr(plhs[0]), root -> rank * root -> rank * sizeof(double));

 mxDestroyArray(plhs[0]);
 mxDestroyArray(plhs[1]);

 return error;
}					/* jacobian */
