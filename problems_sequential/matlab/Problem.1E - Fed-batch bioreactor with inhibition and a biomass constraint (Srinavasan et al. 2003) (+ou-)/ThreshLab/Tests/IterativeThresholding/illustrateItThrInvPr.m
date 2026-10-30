
% This illustrative simulation shows typical behavior for iterative threshold
% methods:
%
%  1) it.dep. ST with argmin(GCVsoft) in each iteration
%  2) it.dep. HT with argmin(GCVhard) in each iteration
%  3) it.dep. ST with argmin(MSEsoft) in each iteration
%  4) it.dep. HT with argmin(MSEhard) in each iteration
%  5) fixed ST with min MSE threshold
%  6) fixed HT with min MSE threshold
%
% typical behaviors:
%  1) large, irregular oscillatory MSE
%  2) smooth, almost monotone decay of MSE
%  3) small, regular oscillatory MSE
%  4) mostly monotone decay of MSE with sudden jumps (change of attraction
%     point in non-convex optimization problem?)
%
% See also
%   simulationItThrInvPr.m
%      which is the same simulation, but without illustrative plots

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
if exist('standardize') ~= 1, standardize = 0; end


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
v0 = vhats0;
msefixedST = zeros(1,maxit);
for it = 1:maxit,
   v1 = iterativeST(w,K,optfixthrS,1,v0);
   msefixedST(it) = norm(v1-v)^2/m;
   v0 = v1;
end
vhatfixedST = v1;

disp('Iteration steps with iteration dependent minimum MSE soft threshold')
v0 = vhats0;
mseitdepSTminmse = zeros(1,maxit);
thritdepSTminmse = zeros(1,maxit);
for it = 1:maxit,
   [v1 thritdepSTminmse(it)] = iterativeMSEST(w,K,v,1,v0,standardize);
   mseitdepSTminmse(it) = norm(v1-v)^2/m;
   v0 = v1;
end
vhatitdepSTminmse = v1;

disp('Iteration steps with iteration dependent minimum GCV soft threshold')
v0 = vhats0;
mseitdepSTmingcv = zeros(1,maxit);
thritdepSTmingcv = zeros(1,maxit);
for it = 1:maxit,
   [v1 thritdepSTmingcv(it)] = iterativeGCVST(w,K,1,v0,standardize);
   mseitdepSTmingcv(it) = norm(v1-v)^2/m;
   v0 = v1;
end
vhatitdepSTmingcv = v1;


figure(1)
plot(msefixedST,'k--','linewidth',3)
hold on
plot(mseitdepSTminmse,'k','linewidth',2)
hold off
title('Iterative Soft Thresholding','fontsize',14,'fontweight','b')
legend('fixed optimal threshold','iteration dependent thr.')
% axis([250,500,0.398,0.404])
a = msefixedST(maxit);
b = min(mseitdepSTminmse);
lo = min(a,b)-abs(b-a)*2;
a = msefixedST(maxit/2);
b = max(mseitdepSTminmse(maxit/2:maxit/2+10));
up = max(a,b) + abs(b-a);
axis([maxit/2,maxit,lo,up])
xlabel('Iteration steps','fontsize',14,'fontweight','b')
ylabel('Loss','fontsize',14,'fontweight','b')
set(gca,'fontsize',14,'fontweight','b')
set(gcf,'position',[0 0 1400 500])
set(gcf, 'PaperPosition', [0 0 14 5]);

figure(2)
plot(mseitdepSTminmse,'k','linewidth',2)
hold on
plot(mseitdepSTmingcv,'r','linewidth',2)
hold off
title('Iterative Soft Thresholding with iteration dependent thresholds',...
      'fontsize',14,'fontweight','b')
legend('minimum MSE thresholds','minimum GCV thresholds')
xlabel('Iteration steps','fontsize',14,'fontweight','b')
ylabel('Loss','fontsize',14,'fontweight','b')
set(gca,'fontsize',14,'fontweight','b')
set(gcf,'position',[0 0 1400 500])
set(gcf, 'PaperPosition', [0 0 14 5]);

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
v0 = vhats0;
msefixedHT = zeros(1,maxit);
for it = 1:maxit,
   v1 = iterativeHT(w,K,optfixthrH,1,v0);
   msefixedHT(it) = norm(v1-v)^2/m;
   v0 = v1;
end
vhatfixedHT = v1;

disp('Iteration steps with iteration dependent minimum MSE hard threshold')
v0 = vhats0;
mseitdepHTminmse = zeros(1,maxit);
thritdepHTminmse = zeros(1,maxit);
for it = 1:maxit,
   [v1 thritdepHTminmse(it)] = iterativeMSEHT(w,K,v,1,v0,standardize);
   mseitdepHTminmse(it) = norm(v1-v)^2/m;
   v0 = v1;
end
vhatitdepHTminmse = v1;

disp('Iteration steps with iteration dependent minimum GCV hard threshold')
v0 = vhats0;
mseitdepHTmingcv = zeros(1,maxit);
thritdepHTmingcv = zeros(1,maxit);
for it = 1:maxit,
   [v1 thritdepHTmingcv(it)] = iterativeGCVHT(w,K,1,v0,standardize);
   mseitdepHTmingcv(it) = norm(v1-v)^2/m;
   v0 = v1;
end
vhatitdepHTmingcv = v1;


figure(3)
plot(msefixedHT,'k--','linewidth',3)
hold on
plot(mseitdepHTminmse,'k','linewidth',2)
hold off
title('Iterative Hard Thresholding','fontsize',14,'fontweight','b')
legend('fixed optimal threshold','iteration dependent thr.')
% axis([250,500,0.398,0.404])
xlabel('Iteration steps','fontsize',14,'fontweight','b')
ylabel('Loss','fontsize',14,'fontweight','b')
set(gca,'fontsize',14,'fontweight','b')
set(gcf,'position',[0 0 1400 500])
set(gcf, 'PaperPosition', [0 0 14 5]);

figure(4)
plot(mseitdepHTminmse,'k','linewidth',2)
hold on
plot(mseitdepHTmingcv,'r','linewidth',2)
hold off
title('Iterative Hard Thresholding with iteration dependent thresholds',...
      'fontsize',14,'fontweight','b')
legend('minimum MSE thresholds','minimum GCV thresholds')
xlabel('Iteration steps','fontsize',14,'fontweight','b')
ylabel('Loss','fontsize',14,'fontweight','b')
set(gca,'fontsize',14,'fontweight','b')
set(gcf,'position',[0 0 1400 500])
set(gcf, 'PaperPosition', [0 0 14 5]);

