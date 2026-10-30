
function f = IWT_BIO(w,L,h,g,sh,sg)

% IWT_BIO - ThreshLab/WT/ -- Inverse fast wavelet transform
%   Usage
%     f = IWT_BIO(w,L,h,g,sh,sg)
%   Inputs
%     w        wavelet coefficients
%     L        coarsest level of transform
%     h,sh     Low pass filter and shift
%     g,sg     High pass filter and shift
%   Outputs
%     f        Inverse wavelet transform Biorthogonal filter
%   Description
%      1D Inverse wavelet transform
%   Note
%     When x is a matrix, an inverse WT is performed on all rows
%
%     sh and sg are shifts taking into account that nonzero segments of
%     filters do not necessarily start at index 1; see makeWTfilter.m and
%     mirrorshifts.m for creation of correct shifts
%
%     see FWT_BIO.m for forward transform
%
%   See also:
%
%     help IWT_BIO
%     help mirrorshifts
%     help makeWTfilter
%


N = size(w,2);
transpose = 0;
if (N == 1)
   N = size(w,1);
   w = w';
   transpose = 1;
end
J = ceil(log2(N));
Ns = N/(2^(J-L));
Nw = Ns;

ph = length(h);
pg = length(g);

%%% SPLITTING THE FILTERS %%%
heven = h(1:2:ph);
hodd  = h(2:2:ph);
geven = g(1:2:pg);
godd  = g(2:2:pg);
   
p = div(max(ph,pg)+1,2)*2;
% p must be even for two reasons:
% 1) it is not, vectors of unequal lengths are added)
% 2) also we are adding a shift of p/2

sho = mod(sh,2);
sgo = mod(sh,2);
sh = div(sh,2);
sg = div(sg,2);
for j = L:J-1,
   s = w(:,1:Ns);
   wj = w(:,Ns+1:Ns+Nw);
   % if-statements in following lines are necessary for matrix inputs and L=0:
   % shift doesn't work as intented in that case. (has undesired effect of
   % shifting in the wrong dimension, since the intended dimension is only one
   % and is therefore not recognized as the dimension in which to perform the
   % shift --- if-statement prevents misinterpretation by shift-command)
   if Ns > 1, s = shift(s,-sh); end
   if Nw > 1, wj = shift(wj,-sg); end
   % periodic extension
   Nj = size(s,2);
   while (2*p > size(s,2))
      s = [s s(:,1:Nj)];
      wj = [wj wj(:,1:Nj)];
   end
   Nj = size(s,2);
   s = [s s(:,1:2*p)];
   wj = [wj wj(:,1:2*p)];
   Njp = size(wj,2);
   x = zeros(size(s,1),2*Njp);
   x1 =  filter(heven,1,s')';
   x2 = filter(geven,1,wj')';
   x3 = filter(hodd,1,s')';
   x4 = filter(godd,1,wj')';
   x(:,1:2:2*Njp) = (1-sho)*x1 + (1-sgo)*x2 + shift(sho*x3 + sgo*x4,-1);
   x(:,2:2:2*Njp) = (1-sho)*x3 + (1-sgo)*x4 + sho*x1 + sgo*x2;
   x(:,1:2*p) = x(:,2*Nj+1:2*Nj+2*p);
   w(:,1:Ns+Nw) = x(:,1:Ns+Nw);
   Ns = Ns+Nw;
   Nw = Ns;
end
f = w;
if transpose, f = f'; end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
