
function lambda = RCVS(yp,M,h,ht)

% usage: lambda = RCVS(yp,M,h,ht);
%
% Redundant Conditional Variance Stabilisation for Poisson intensity estimation
%
% Non-decimated wavelet transform, Cross Validation
%
% output: lambda: vector of intensities
% input : yp    : Poisson counts, i.e., the observations
%         L,h,ht: wavelet transform parameters, see RWT_BIO.m, IRT_BIO.m
%

N = length(yp);
if nargin < 4,
   ht = h;
end
[h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);
w = RWT_BIO(yp,M,ht,gt,sht,sgt);
pn  = RWTpoisnorm(yp,M,ht,gt,sht,sgt);
pn = pn + (pn==0);
wpn = w./pn;

wpnthr = wpn;
for j = 1:M,
   j
   % select coefficients at scale j
   ww = wpn(j,1:N);
   [wwthr thr] = mingcv(ww);
   wpnthr(j,:) = wwthr;
end
wthr = wpnthr.*pn;
lambda = IRT_BIO(wthr,h,g,sh,sg);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
