
% This test compares CVS with Anscombe's method on a single simulation.
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
  testsignal = bumps(x); 
  testsignal = 10*(testsignal+1);
end
f = testsignal;
yp = poisson(f,realrandom);

if ~exist('h'), % wavelet filters
   [h, ht] = makeWTfilter('Daubechies',3);
   [h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);
end
if ~(exist('L') == 1), L = 5; end % lowest resolution level
J = ceil(log2(N));
M = J-L;

% --1-- CVS - wavelet transform

w = RWT_BIO(yp,M,ht,gt,sht,sgt);
pn = RWTpoisnorm(yp,M,ht,gt,sht,sgt);
pn = pn + (pn==0);
wpn = w./pn;
v = RWT_BIO(f,M,ht,gt,sht,sgt);
vpn = v./pn;

vvhat = wpn;
for j = 1:M,
   ww = wpn(j,1:N);
   vv = vpn(j,1:N);
   [wwthr thr] = mingcv(ww);
   [wwthr thr] = minmse(vv,ww);
   vvhat(j,1:N) = wwthr;
end
vhat = vvhat.*pn;
fhat = IRT_BIO(vhat,h,g,sh,sg);

% --2-- Anscombe...

if ~exist('c'),
   c = 3/8;
end
wans = RWT_BIO(sqrt(yp+c),M,ht,gt,sht,sgt);
vans = RWT_BIO(sqrt(f+c),M,ht,gt,sht,sgt);
vhatans = wans;
for j = 1:M,
   ww = wans(j,1:N);
   vv = vans(j,1:N);
   [wwthr thr] = mingcv(ww);
   [wwthr thr] = minmse(vv,ww);
   vhatans(j,1:N) = wwthr;
end
fhatans = IRT_BIO(vhatans,h,g,sh,sg);
fhatans = fhatans.^2-c;
fhatans = fhatans + VsqrtPoisson(fhatans,c);
% VsqrtPoisson(fhatans) is an estimation of the additional bias due to
% sqrt+square operations

snrin = snr(f,yp-f);
snrcvs = snr(f,fhat-f);
snrans = snr(f,fhatans-f);

L1errcvs = mean(abs(fhat-f));
L1errans = mean(abs(fhatans-f));
L2errcvs = mean((fhat-f).^2);
L2errans = mean((fhatans-f).^2);

biascvs = mean(fhat-f);
biasans = mean(fhatans-f);

