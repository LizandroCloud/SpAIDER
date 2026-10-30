
function wcoef = RWT_BIO(x,M,ht,gt,sht,sgt)

% RWT_BIO - ThreshLab/WT/ -- Redundant 1D Biorthogonal WT
%  Usage
%    w = RWT_BIO(x,M,ht,gt,sht,sgt)
%  Inputs
%    x         signal
%    M         number of levels
%    ht        Low pass analysis filter
%    gt        High pass analysis filter
%    sht       Low pass analysis filter shift (see below)
%    sgt       High pass analysis filter shift (see below)
%  Outputs
%    w         the wavelet coefficients; (M+1) x length(x) matrix
%              j-th row has wavelet/detail coefficients at level j
%              (M+1)st row has coarse scale (=level J-M) scaling coef
%  Description
%    Computes Redundant 1D Biorthogonal WT
%    The redundant wavelet transform is known under many names including
%       translation-invariant wavelet transform
%       stationary wavelet transform
%       a-trous algorithm (although this term is ambivalent)
%       cycle spinning
%       non-decimated wavelet transform
%  Note
%    sht and sgt are shifts taking into account that nonzero segments of
%    filters do not necessarily start at index 1; see makeWTfilter.m and
%    mirrorshifts.m for creation of correct shifts
%
%    see IRT_BIO.m for inverse transform
%
%  See also:
%    help FWT_BIO
%    help IRT_BIO (inverse redundant transform)
%    help mirrorshifts
%    help makeWTfilter


N = length(x);
Nj = N;
pht = length(ht);
pgt = length(gt);
sht = -sht - pht + 1;
sgt = -sgt - pgt + 1;
   
% filters must be flipped in forward wavelet transform
% flipping causes change in shifts sht and sgt (see above)
ht = ht(pht:-1:1);
gt = gt(pgt:-1:1);
p = div(max(pht,pgt)+1,2)*2;
% p must be even for two reasons:
% otherwise vectors of dissimilar lengths are being added
% and also we are adding a shift of p/2

w = zeros(M+1,N);
% w(M+1,:,:) = x;
w(M+1,:) = row(x);
% last row of w will contain scaling coefficients (on coarsest scale)
% x is being recycled as working variable

skip = 1;
for j = 1:M
  s = shift(w(M+1,:),-skip*sgt);
  % periodic exetension
  while (p*skip > size(s,2))
      for tel = 1:skip
         s = [s s(:,1:N)];
      end
  Nj = size(s,2);
  end
  s = [s s(:,1:skip*p)];
  wj = zeros(size(s));
  Njp = size(wj,2);
  for k = 1:skip
      sk = s(:,k:skip:Njp);
      wj(k:skip:Njp) = filter(gt,1,sk);
  end
  s = shift(w(M+1,:),-skip*sht);
  % periodic exetension
  while (p*skip > size(s,2))
      for tel = 1:skip
         s = [s s(:,1:N)];
      end
  Nj = size(s,2);
  end
  s = [s s(:,1:skip*p)];
  x = zeros(size(s));
  Njp = size(x,2);
  for k = 1:skip
      sk = s(:,k:skip:Njp);
      x(k:skip:Njp) = filter(ht,1,sk);
  end
  wj(:,1:skip*p) = wj(:,Nj+1:Nj+skip*p);
  x(:,1:skip*p) = x(:,Nj+1:Nj+skip*p);
  w(j,:) = wj(:,1:N);
  w(M+1,:) = x(:,1:N);
  skip = skip*2;
end
wcoef = w;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
