
function w = shiftCWT(w,M,m);

% shiftCWT - ThreshLab2/ContinuousUnbalancedHaar/
%   Usage
%     w = shiftCWT(w,M,m);
%   Inputs
%     w   continuous wavelet transform coefficients
%     M   number of dyadic levels/scales (depth of transform)
%     m   number of intermediate scales (see help CHT)
%   Outputs
%     w   shifted continuous wavelet transform coefficients
%   Description
%     shifts wavelet coefficients such that transform has centralised basis
%     functions (instead of left allined basis functions)
%   See also
%      help CHT


Mm = size(w,1);

r = 0;
mm = 0;
while mm < m,
   r = r+1;
   mm = 2^(r+1)-r-2;
end
r = min(r,M-2);
m = 2^(r+1)-r-2;

s = 1/2;
for j = 2:M-r, % dyadic shifts
   s = s*2;
   w(j,:) = shift(w(j,:),-s);
end
d = s;
for j = M-r+1:Mm, % linear shifts
   s = s+d;
   w(j,:) = shift(w(j,:),-s);
end
