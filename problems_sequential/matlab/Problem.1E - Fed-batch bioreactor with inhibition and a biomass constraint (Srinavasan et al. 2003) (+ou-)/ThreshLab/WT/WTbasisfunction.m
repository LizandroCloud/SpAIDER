
function psi = WTbasisfunction(M,h,g,sh,sg,scaling)

% usage: psi = WTbasisfunction(M,h,g,sh,sg,scaling);
%  or  : psi = WTbasisfunction(M,ht,gt,sht,sgt,scaling):
% (i.e., procedure possible for inverse and analysis basis functions)
%
% M is number of resolution levels. M = J-L;
% scaling is optional boolean (default 0). If it is set to 1, this function
% computes scaling basis, otherwise wavelet basis.

if nargin < 6,
   scaling = 0;
end
nh = length(h);
ng = length(g);
rjs = 1; 
for j = 1:M,
   rjw = rjs+2^(j-1)*ng-2^(j-1)+1-1;
   rjs = rjs+2^(j-1)*nh-2^(j-1)+1-1;
end
if scaling,
   J = ceil(log2(rjs));
   N = 2^J;
   L = J-M;
   i = ceil(2^(L+1)/4);
else
   J = ceil(log2(rjw));
   N = 2^J;
   L = J-M;
   i = ceil(3*2^(L+1)/4);
end
x0 = zeros(1,N);
x0(i) = 1;
gamma = IWT_BIO(x0,L,h,g,sh,sg);
x = zeros(1,N);
[first, last] = FWTdualsupport(x,L,h,g,sh,sg);
if first(i) > last(i),
   gamma = [gamma(first(i):N) gamma(1:last(i))];
else
   gamma = gamma(first(i):last(i));
end
psi = gamma;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
