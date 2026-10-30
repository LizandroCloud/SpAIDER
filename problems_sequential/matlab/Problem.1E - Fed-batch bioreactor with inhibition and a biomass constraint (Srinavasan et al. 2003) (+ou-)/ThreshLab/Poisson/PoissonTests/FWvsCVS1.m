
% Fisz Wavelet (Fryzlewicz, e.a., 2004)
% versus Conditional Variance Stabilisation

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
   [h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);
end
if ~(exist('L') == 1), L = 5; end % lowest resolution level
J = ceil(log2(N));

% --1-- CVS

w = FWT_BIO(yp,L,ht,gt,sht,sgt);
pn = FWTpoisnorm(yp,L,ht,gt,sht,sgt);
pn = pn + (pn==0);
wpn = w./pn;
v = FWT_BIO(f,L,ht,gt,sht,sgt);
vpn = v./pn;

vvhat = wpn;
for j = L:J-1,
   ww = wpn(2^j+1:2^(j+1));
   vv = vpn(2^j+1:2^(j+1));
   [wwthr thr] = mingcv(ww);
   [wwthr thr] = minmse(vv,ww);
   vvhat(2^j+1:2^(j+1)) = wwthr;
end
vhat = vvhat.*pn;
fhat = IWT_BIO(vhat,L,h,g,sh,sg);

% --2-- Haar-Fisz-Wavelet

whfw = FFHFW(yp,L,ht,gt,sht,sgt);
vhfw = FFHFW(f,L,ht,gt,sht,sgt);
vhathfw = whfw;
for j = L:J-1,
   ww = whfw(2^j+1:2^(j+1));
   vv = vhfw(2^j+1:2^(j+1));
   [wwthr thr] = mingcv(ww);
   [wwthr thr] = minmse(vv,ww);
   vhathfw(2^j+1:2^(j+1)) = wwthr;
end
fhathfw = IFHFW(vhathfw,L,h,g,sh,sg);

% comparison

snrin = snr(f,yp-f);
snrcvs = snr(f,fhat-f);
snrhfw = snr(f,fhathfw-f);

L1errcvs = mean(abs(fhat-f));
L1errhfw = mean(abs(fhathfw-f));
L2errcvs = mean((fhat-f).^2);
L2errhfw = mean((fhathfw-f).^2);

biascvs = mean(fhat-f);
biashfw = mean(fhathfw-f);

