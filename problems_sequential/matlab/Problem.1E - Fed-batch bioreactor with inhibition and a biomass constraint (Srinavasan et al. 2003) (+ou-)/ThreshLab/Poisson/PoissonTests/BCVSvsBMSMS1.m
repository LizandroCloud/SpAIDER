
% This test compares BCVS with BMSMS (Kolaczyk, 1999) on a single simulation.
% non-decimated wavelet transform

if ~exist('realrandom'), realrandom = 0; end % binary value indicating whether
                                             % or not the noise generator
                                             % starts from a fixed seed.
if ~exist('x'),
   if ~exist('N'), N = 2048; end
   x = makegrid(N,0);
else
   N = length(x);
end
if ~exist('testsignal'),
  testsignal = heavisine(x);
  testsignal = testsignal - min(testsignal);
  testsignal = 5*testsignal+1;
end
f = testsignal;
yp = poisson(f,realrandom);

if ~exist('h'), % wavelet filters
   [h, ht] = makeWTfilter('Daubechies',3);
end
if ~(exist('L') == 1), L = 5; end % lowest resolution level
J = ceil(log2(N));
M = J-L;

% --1-- BRCVS Bayesian Redundant Conditional Variance Stabilisation

yshr = BRCVS(yp,M,h,ht);

% --2-- BMSMS, Kolaczyk JASA, 1999, pp 920 - 933

if ~exist('empirical'),
   empirical = 0;
end
if exist('BMSMShrink') ~= 2,
   disp('--------------------------------------------------------------')
   disp('WARNING: Cannot find BMSMShrink.m')
   disp('Make sure to include the routine BMSMShrink.m into your path.')
   disp('This function is part of the software made available by ')
   disp('T. Sapatinas and collaborators. See')
   disp('http://www.ucy.ac.cy/~fanis/links/software.html');
   disp('--------------------------------------------------------------')
   yshrBMSMS = yp;
else
   [yshrBMSMS pBMSMS aBMSMS] = BMSMShrink(yp,empirical+1);
end


snrin = snr(f,yp-f);
snrbcvs = snr(f,yshr-f);
snrbmsms = snr(f,yshrBMSMS-f);

L1errbcvs = mean(abs(yshr-f));
L1errbmsms = mean(abs(yshrBMSMS-f));
L2errbcvs = mean((yshr-f).^2);
L2errbmsms = mean((yshrBMSMS-f).^2);

biasbcvs = mean(yshr-f);
biasbmsms = mean(yshrBMSMS-f);

