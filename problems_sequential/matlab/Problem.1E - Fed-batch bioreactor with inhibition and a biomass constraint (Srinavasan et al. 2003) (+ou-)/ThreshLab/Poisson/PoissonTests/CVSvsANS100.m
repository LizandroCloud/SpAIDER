
% This file runs a number of (100 by default) simulations of CVSvsANS1.m
% see also CVSvsANS1.m

realrandom = 1; % binary value indicating whether or not the noise generator
                % starts from a fixed seed.
rand('state',12345);
randn('state',12345);
% initial seeds

if ~exist('x'),
   if ~exist('N'), N = 2048; end
   x = makegrid(N,0);
else
   N = length(x);
end
if ~exist('testsignal'),
  testsignal = bumps(x);
  testsignal = 10*(testsignal+1);
end

if ~exist('h'), % wavelet filters
   [h, ht] = makeWTfilter('Daubechies',3);
   [h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);
end
if ~(exist('L') == 1), L = 5; end % lowest resolution level

if ~exist('samplesize'), samplesize = 100; end

for iter = 1:samplesize,
   iter
   CVSvsANS1
   SNRans1(iter) = snrans;
   SNRcvs1(iter) = snrcvs;
end

% flip testsignal
testsignal = max(testsignal) + min(testsignal) - testsignal;

for iter2 = 1:samplesize,
   iter2
   CVSvsANS1
   SNRans2(iter2) = snrans;
   SNRcvs2(iter2) = snrcvs;
end

disp('table values:')
meanSNRcvs1 = mean(SNRcvs1)
meanSNRans1 = mean(SNRans1)
meanSNRcvs2 = mean(SNRans2)
meanSNRans2 = mean(SNRcvs2)

figure(1)
plot(x,testsignal)
title('Intensity function')
figure(2)
plot(x,yp)
title('Typical realisation')
