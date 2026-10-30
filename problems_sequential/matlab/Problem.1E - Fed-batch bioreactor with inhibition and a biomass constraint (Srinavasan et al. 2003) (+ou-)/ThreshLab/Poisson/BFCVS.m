
function lambda = BFCVS(yp,L,h,ht)

% usage: lambda = BFCVS(yp,L,h,ht);
%
% Bayesian Fast (= decimated) Conditional Variance Stabilisation for 
% Poisson intensity estimation
%
% output: lambda: vector of intensities
% input : yp    : Poisson counts, i.e., the observations
%         L,h,ht: wavelet transform parameters, see RWT_BIO.m, IRT_BIO.m
%
% This method is empirical Bayes: it includes a maximum marginal likelihood 
% estimation of mixture parameter p, see MmLEmixnorm.m
%
% The model requires simultanuous access to both wavelet coefficients and input
% observations, so there is no elegant implementation of the form:
% w = WT(input)
% wshr = shrink(w)
% output = IWT(wshr)
% The second step is not based on wavelet coefficients only, this is why we
% decide to bring the whole procedure into a single function, this one.

N = length(yp);
J = log2(N);
if nargin < 4,
   ht = h;
end
[h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);
w = FWT_BIO(yp,L,ht,gt,sht,sgt);
[pn first last] = FWTpoisnorm(yp,L,ht,gt,sht,sgt);
% the vectors first and last contain for each wavelet coefficient
% the index of the first and last element of the input vector involved in the 
% computation of that coefficient
pn = pn + (pn==0);
wpn = w./pn;

a = 1; % parameter of dirichlet model, could also be estimated empirically (see
       % paper)

wpnshr = wpn;
for j = L:J-1,
   j
   % select coefficients at scale j
   ww = wpn(2^j+1:2^(j+1));
   pnj = pn(2^j+1:2^(j+1));
   Nvec = pnj.^2; % poisson normalisation pn equals sqrt(Nvec), Nvec is sum of
                  % observations involved in computation of wavelet coefficient
   wwshr = ww;
   
   firstj = first(2^j+1:2^(j+1));
   lastj = last(2^j+1:2^(j+1));
   rj = lastj-firstj+1;
   rj = rj + N*(rj<0);

   gammaj = WTbasisfunction(J-j,ht,gt,sht,sgt);
   sig0j = sqrt(norm(gammaj)^2/rj(1));
   sig0jest = MAD1(ww)/0.6745;
   if J-j<7, sig0j = sig0jest; end

   A = a*rj;
   sig1j = sig0j*sqrt((Nvec+A)./(1+A)); % note that sig1j is a
                                             % non-constant vector: the
                                             % standard deviation of
                                             % significant coefficients is not
                                             % a constant
   % V(ww) = p(j) * sig1j^2 + (1-p(j)) * sig0j^2

   p(j) = MmLEmixnorm(ww,sig0j,sig1j);
   n = length(ww);
   if j > 1, p(j) = max(p(j),p(j-1)); end
   C = sqrt((Nvec+A)./(1+A));
   D = (Nvec-1)./(Nvec+A);
   pmin(2^j+1:2^(j+1)) = C./(C+n.^D);
   pprior = max(pmin(2^j+1:2^(j+1)),p(j));
   pvec(2^j+1:2^(j+1)) = pprior;
   qprior = 1 - pprior;
   ppost = pprior.*gauss(ww,sig1j)./...
           (pprior.*gauss(ww,sig1j) + qprior.*gauss(ww,sig0j));
   mu1postj = ww.*(Nvec./(Nvec+A));

   % mu1postj is posterior mean, now proceed to posterior median...
   % The posterior mean depends on input value and normalisation only. The
   % posterior distribution (including the variance) as a whole however depends
   % on all observations involved in the computation of that coefficient (see
   % paper)
 
   for i = 1:n
      if firstj(i) > lastj(i),
         xx = [yp(firstj(i):N) yp(1:lastj(i))];
      else
         xx = yp(firstj(i):lastj(i));
      end
      nn = sum(xx);
      % nn = Nvec(i);
      aa = (a+xx)/(A(i)+nn);
      sig1post = sum(gammaj.^2.*aa)-(sum(gammaj.*aa))^2;
      sig1post = sig1post*nn/(nn+A(i)+1);
      sig1postj(i) = sig1post;
      wwshr(i) = medianmass0norm(ppost(i),mu1postj(i),sig1post);
   end
   wpnshr(2^j+1:2^(j+1)) = wwshr;
end
wshr = wpnshr.*pn;
lambda = IWT_BIO(wshr,L,h,g,sh,sg);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
