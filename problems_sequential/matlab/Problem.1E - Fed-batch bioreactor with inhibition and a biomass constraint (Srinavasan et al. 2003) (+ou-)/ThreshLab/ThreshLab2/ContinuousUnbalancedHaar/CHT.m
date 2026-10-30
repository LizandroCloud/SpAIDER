
function [w, v] = CHT(y,M,m);

% CHT - ThreshLab2/ContinuousUnbalancedHaar/ -- Continuous Haar Transform
%   Usage
%     [w, v] = CHT(y,M,m);
%   Inputs
%     y   observed values of function to be analysed
%     M   number of dyadic levels/scales (depth of transform)
%     m   number of intermediate scales (see below; default is 0)
%   Outputs
%     w   continuous wavelet transform coefficients
%     v   scaling coefficients at each scale
%   Description
%     evaluates the Continuous Haar Transform in a finite number of scales,
%     specified by parameters M and m. M is the number of dyadic scales, as in
%     the Discrete Wavelet Transform (decimated or undecimated). 
%     m is number of additional, intermediate scales; m is ceiled to the next
%     value of 2^(r+1)-r-2 for integer r. This is to make an equidistant
%     distribution of the additional scales in between the coarsest dyadic
%     scales.
%     If m = 0, then this function coincides with the undecimated Haar wavelet
%     transform, see  RWT_BIO
%   See also
%      help RWT_BIO
%      help FWT_BIO

if nargin < 3,
   m = 0;
end

N = length(y);
s = zeros(M,N);
s(1,:) = row(y);
for j = 2:M,
   s(j,:) = s(j-1,:) + shift(s(j-1,:),2^(j-2));
end

r = 0;
mm = 0;
while mm < m,
   r = r+1;
   mm = 2^(r+1)-r-2;
end
r = min(r,M-2);
m = 2^(r+1)-r-2;

Mm = M+m;
inbetween = [ones(1,M-r) 2.^(1:r)];
dyadicscales = cumsum(inbetween);

for j = 1:M,
   w(dyadicscales(j),:) = s(j,:) - shift(s(j,:),2^(j-1));
   v(dyadicscales(j),:) = s(j,:) + shift(s(j,:),2^(j-1));
end
for j = M-r:M-1,
   s0 = s(j,:);
   s1 = shift(s(M-r-1,:),2^(j-1));
   for jj = 1:2^(j-M+r+1)-1,
      s0 = s0+s1;
      w(dyadicscales(j)+jj,:) = s0 - shift(s0,2^(j-1)+jj*2^(M-r-2));
      v(dyadicscales(j)+jj,:) = s0 + shift(s0,2^(j-1)+jj*2^(M-r-2));
      s1 = shift(s1,2^(M-r-2));
   end
end
