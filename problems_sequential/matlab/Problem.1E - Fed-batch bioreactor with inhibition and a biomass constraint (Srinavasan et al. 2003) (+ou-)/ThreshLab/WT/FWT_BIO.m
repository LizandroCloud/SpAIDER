
function wcoef = FWT_BIO(x,L,ht,gt,sht,sgt)

% FWT_BIO - ThreshLab/WT/ -- Fast 1D Biorthogonal WT
%  Usage
%    wcoef = FWT_BIO(x,L,ht,gt,sht,sgt)
%  Inputs
%    x         signal
%    L         coarsest level
%    ht        Low pass analysis filter
%    gt        High pass analysis filter
%    sht       Low pass analysis filter shift (see below)
%    sgt       High pass analysis filter shift (see below)
%  Outputs
%    wcoef     the wavelet coefficients
%  Description
%    Computes Fast 1D Biorthogonal WT
%  Note
%    When x is a matrix, a WT is performed on all rows
%
%    sht and sgt are shifts taking into account that nonzero segments of 
%    filters do not necessarily start at index 1; see makeWTfilter.m and
%    mirrorshifts.m for creation of correct shifts
%
%    see IWT_BIO.m for inverse transform
%
%  See also:
%
% help IWT_BIO
% help mirrorshifts
% help makeWTfilter
%

pht = length(ht);
pgt = length(gt);
sht = sht + pht - 1;
sgt = sgt + pgt - 1;
   
% filters must be flipped in forward wavelet transform
ht = ht(pht:-1:1);
gt = gt(pgt:-1:1);
p = div(max(pht,pgt)+1,2)*2;
% p must be even for two reasons:
% 1) it is not, vectors of unequal lengths are added)
% 2) also we are adding a shift of p/2

%%% SPLITTING THE FILTERS %%%
hteven = ht(1:2:pht);
htodd  = ht(2:2:pht);
gteven = gt(1:2:pgt);
gtodd  = gt(2:2:pgt);

N = size(x,2);
transpose = 0;
if (N == 1)
   x = row(x);
   N = size(x,2);
   transpose = 1;
end
w = x;
s = x;
Ns = N;

J = ceil(log2(N));
for j=J-1:-1:L
   Nj = size(s,2);
   Ns = div(Nj+1,2);
   Nw = div(Nj,2);
   % periodic extension
   while (p > size(s,2))
      s = [s s(:,1:Nj)];
   end
   Nj = size(s,2);
   s = shift(s,sgt);
   s = [s s(:,1:p)];

   Njp = size(s,2);
   Nj = div(Nj+1,2);
 
   %%% SPLITTING THE SIGNAL (=> downsampling)
   seven = s(:,1:2:Njp);
   sodd  = s(:,2:2:Njp);
   wj = filter(gteven,1,seven')' + shift(filter(gtodd,1,sodd')',-1);
   s = s(:,1:Njp-p);
   s = shift(s,sht-sgt);
   s = [s s(:,1:p)];
   seven = s(:,1:2:Njp);
   sodd  = s(:,2:2:Njp);
   s = filter(hteven,1,seven')' + shift(filter(htodd,1,sodd')',-1);
  
   % matlab filter.m assumes causal signals. Since we work with periodic
   % extensions, the first coefficients are wrong. Therefore, we need a shift of
   % the input (dropping the first coefficients after shifting)
   wj(:,1:p/2) = wj(:,Nj+1:Nj+p/2);
   s(:,1:p/2) = s(:,Nj+1:Nj+p/2);
   wj = wj(:,1:Nw);
   s = s(:,1:Ns);
   w(:,Ns+1:Ns+Nw) = wj;
   clear wj
end
w(:,1:Ns) = s(:,1:Ns);
wcoef = w;
if transpose, wcoef = wcoef'; end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
