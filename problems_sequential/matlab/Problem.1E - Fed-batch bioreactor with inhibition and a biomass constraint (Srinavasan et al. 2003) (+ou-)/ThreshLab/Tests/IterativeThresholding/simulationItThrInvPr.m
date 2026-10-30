
% simulation of iterative thresholding routines
%  See also illustrateItThrInvPr.m, which runs the same simulations, with
%  illutsrative plots
%  See also simulstudyItThrInvPr.m , which is a study of several runs of this
%  file

if exist('typeofK') ~= 1, typeofK = 'random'; end
if exist('realrandom') ~= 1, realrandom = 0; end
if exist('dimK') ~= 1, dimK = 0; end
   if ~dimK, m = 1000; n = 200; end % n = number of observations
clear dimK
if exist('snratio') ~= 1, snratio = 10; end
if exist('degreeofsparsity') ~= 1, degreeofsparsity = 0.05; end
if exist('bandwidth') ~= 1, bandwidth = 1; end
if exist('singvalues') ~= 1, singvalues = 1; end
if exist('n_singval') ~= 1, n_singval = 1; end


p = degreeofsparsity;

[w K v stdev Z] = ...
  setupsimulationYisKxplussigmaZ(typeofK,m,n,p,snratio,realrandom,...
                                'bandwidth',bandwidth,...
                                'singularvalues',singvalues,...
                                'n_singval',n_singval);
Kv = K*v;
[n m] = size(K);

%%%%%%%%%%%%%%%% SOFT THRESHOLDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('1. ITERATIVE SOFT THRESHOLING')
% first find minimum MSE fixed threshold
warning off all

initit = 100;
thr = sqrt(2*log(m))*stdev/3;
vhats0 = iterativeST(w,K,thr,initit);

maxit = 500;
tmin = 0;
tmax = sqrt(2*log(m))*sqrt(var2(column(w)-K*column(v)));
optfixthrS = minfunc1Dfibon('MSEitthresh','argpos',4,'lo',tmin,'up',tmax,...
                             v,w,K,maxit,vhats0,1);
% [mse vhat] = MSEitthresh(v,w,K,t,maxit,vhat0,softORhard);

disp('Iteration steps with fixed soft threshold')
v1fixedST = iterativeST(w,K,optfixthrS,maxit,vhats0);

disp('Iteration steps with iteration dependent minimum MSE soft threshold')
v1itdepSTminmse = iterativeMSEST(w,K,v,maxit,vhats0);

disp('Iteration steps with iteration dependent minimum GCV soft threshold')
v1itdepSTmingcv = iterativeGCVST(w,K,maxit,vhats0);

%%%%%%%%%%%%%%%% HARD THRESHOLDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('2. ITERATIVE HARD THRESHOLING')

% first find minimum MSE fixed threshold
warning off all

initit = 100;
thr = sqrt(2*log(m))*stdev/3;
vhats0 = iterativeST(w,K,thr,initit);

maxit = 500;
tmin = 0;
tmax = sqrt(2*log(m))*sqrt(var2(column(w)-K*column(v)));
optfixthrH = minfunc1Dfibon('MSEitthresh','argpos',4,'lo',tmin,'up',tmax,...
                             v,w,K,maxit,vhats0,0);
% [mse vhat] = MSEitthresh(v,w,K,t,maxit,vhat0,softORhard);

disp('Iteration steps with fixed hard threshold')
v1fixedHT = iterativeHT(w,K,optfixthrH,maxit,vhats0);

disp('Iteration steps with iteration dependent minimum MSE hard threshold')
v1itdepHTminmse = iterativeMSEHT(w,K,v,maxit,vhats0);

disp('Iteration steps with iteration dependent minimum GCV hard threshold')
v1itdepHTmingcv = iterativeGCVHT(w,K,maxit,vhats0);
