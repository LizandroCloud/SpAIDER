
% This test illustrates the application of GCV as a variable selection
% criterion in the model Y = K*beta + stdev*Z, with Z ~ i.i.d. N(0,1)
% the variable selection proceeds through optimization of an ell_1 penalized
% least squares (LASSO), implemented by iterative soft thresholding.
%
% Several design matrices are possible, see
%      help setupsimulationYisKxplussigmaZ
%

if exist('typeofK') ~= 1, typeofK = 'bandrandom'; end
if exist('realrandom') ~= 1, realrandom = 0; end
if exist('dimKgiven') ~= 1, dimKgiven = 0; end
   if ~dimKgiven, m = 3000; n = 1000; end % m = number of variables
                                          % n = number of observations
clear dimKgiven
if exist('snratio') ~= 1, snratio = 10; end
if exist('degreeofsparsity') ~= 1, degreeofsparsity = 0.01; end
if exist('bandwidth') ~= 1, bandwidth = 20; end
if exist('singvalues') ~= 1, singvalues = 1; end
if exist('n_singval') ~= 1, n_singval = min(m,n); end
if exist('standardize') ~= 1, standardize = 0; end
if exist('xvdm') ~= 1, xvdm = (1:n)/n; end

% other options include:
%   typeofK = 'SVDgivenSV';
%   singvalues = (1:m).^2/m.^2;
%
%   typeofK = 'Vandermonde';
%
%   typeofK = 'eye';

p = degreeofsparsity;
if ~realrandom,
   randn('state',3); rand('state',3);
   disp('simulation setup with real random, but state has been fixed before.')
end

[Y K beta stdev Z] = ...
  setupsimulationYisKxplussigmaZ(typeofK,m,n,p,snratio,1,...
                                'bandwidth',bandwidth,...
                                'singularvalues',singvalues,...
                                'n_singval',n_singval,...
                                'vdmcovar',xvdm);
Kbeta = K*beta;
[n m] = size(K);
condKKT = cond(K*K')

if exist('maxit') ~= 1, maxit = 3000; end
if exist('HT')~= 1, HT = 0; end
% HT is a binary variable: HT=1 activates a least squares postprocessing within
% the selected model. (If K = I, this is hard thresholding)
% 

if exist('nthr') ~= 1, nthr = 50; end

% iterative soft-thresholding
if nthr > 0,

univthr = sqrt(2*log(n))*stdev;
if exist('maxthrfactor') ~= 1, maxthrfactor = 1; end
maxthr = univthr*maxthrfactor;
thr = (1:nthr)/nthr*maxthr;
prederrorITST = zeros(1,nthr);
N1 = zeros(1,nthr);
SSeITST = zeros(1,nthr);
for t = 1:nthr,
t
   betahats = iterativeST(Y,K,thr(t),maxit);
   if HT==1, % post-processing: orthogonal projection (least squares)
      X = find(abs(betahats)>eps); Kx = K(1:n,X);
      betahats(X) = (Kx'*Kx)\(Kx'*column(Y));
   end
   prederrorITST(t) = norm(K*betahats-Kbeta)^2/n;
   SSeITST(t) = norm(K*betahats-Y)^2;
   N1(t) = sum(abs(betahats)>eps);
end
CpITSTstdev = SSeITST/n +2*N1/n*stdev^2-stdev^2;
GCVITST = (SSeITST/n ./ (1-N1/n).^2);
[mingcv argmingcv] = min(GCVITST);
[minCpstdev argminCpstdev] = min(CpITSTstdev);
stdevestGCV = sqrt(SSeITST(argmingcv)/(n-N1(argmingcv)));
stdevestCp = sqrt(SSeITST(argminCpstdev)/(n-N1(argminCpstdev)));
% we use the variance estimator of GCV to perform a practical Cp procedure.
stdevest = stdevestGCV;
CpITSTstdevest = SSeITST/n +2*N1/n*stdevest^2-stdevest^2;

plot(thr,prederrorITST,'r','linewidth',2)
hold on 
plot(thr,CpITSTstdev,'b','linewidth',2)
plot(thr,GCVITST-stdev^2,'k','linewidth',2)
hold off
legend('pred error','Cp -- SURE','GCV-\sigma^2')
set(gca,'fontsize',14,'fontweight','bold')

relerr = ((2*N1/n)-(N1/n).^2)./(1-N1/n).^2 - ...
          (N1/n).^2./(1-N1/n).^2.*stdev^2./CpITSTstdev;
abserr = relerr.*CpITSTstdev;

% quality of Cp with exact stdev and with estimated stdev
[minCpstdevest argminCpstdevest] = min(CpITSTstdevest);
dGCVITST = GCVITST(2:nthr)-GCVITST(1:nthr-1);
% instead of the absolute minimum of the GCV curve, we may consider the largest
% local minimizer. This way, the procedure stays away from the region where GCV
% fluctuates due to instabilities near lambda=0.
argmingcv2 = max(find(dGCVITST<-eps))+1;
minPE = min(prederrorITST);
EffGCV        = minPE/prederrorITST(argmingcv)
EffGCV2       = minPE/prederrorITST(argmingcv2)
EffCpstdev    = minPE/prederrorITST(argminCpstdev)
EffCpstdevest = minPE/prederrorITST(argminCpstdevest)
% hold on
% plot(thr,CpITST+abserr,'r--')
% hold off


end % if nthr > 0
