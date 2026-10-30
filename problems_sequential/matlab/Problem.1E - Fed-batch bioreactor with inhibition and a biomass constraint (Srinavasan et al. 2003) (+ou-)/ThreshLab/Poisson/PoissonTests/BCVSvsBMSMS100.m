
% This file runs a number of (100 by default) simulations of BCVSvsBMSMS1.m
% see also BCVSvsBMSMS1.m

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
  testsignal = heavisine(x);
  testsignal = testsignal - min(testsignal);
  testsignal = 5*testsignal+1;
end

if ~exist('h'), % wavelet filters
   [h, ht] = makeWTfilter('Daubechies',3);
end
if ~(exist('L') == 1), L = 5; end % lowest resolution level

if ~exist('empirical'),
   empirical = 0; % if empirical == 1, BMSMS includes empirical estimation of
                  % hyperparameters (EM-algorithm)
end

if ~exist('samplesize'), samplesize = 100; end

for iter = 1:samplesize,
   iter
   BCVSvsBMSMS1
   SNRbmsms1(iter) = snrbmsms;
   SNRbcvs1(iter) = snrbcvs;
end

disp('table values:')
meanSNRbcvs = mean(SNRbcvs1)
meanSNRbmsms = mean(SNRbmsms1)

figure(1)
plot(x,testsignal)
title('Intensity function')
figure(2)
plot(x,yp)
title('Typical realisation')

