
function [first, last] = FWTdualsupport(x,L,ht,gt,sht,sgt)

% usage: [first, last] = FWTdualsupport(x,L,ht,gt,sht,sgt);
%
% This function computes the support of the (discrete) dual wavelet functions,
% i.e. the first and last non-zero entry of the forward wavelet transform
% (output first and last)
% The parameters L,ht,gt,sht,sgt correspond to those in the function call 
% w = FWT_BIO(x,L,ht,gt,sht,sgt);
%
N = length(x);
nht = length(ht);
ngt = length(gt);
first = (1:N);
last = first;
J = ceil(log2(N));
rjs = 1; % rj = last(i)-first(i) (if first(i) < last(i), otherwise modulo N)
for j = J-1:-1:L,
   rjw = rjs+2^(J-1-j)*ngt-2^(J-1-j)+1-1;
   rjs = rjs+2^(J-1-j)*nht-2^(J-1-j)+1-1;
   % reduces to: rjs = 2*rjs+nht-2
   % to understand the previous two statements, one has to realise that rjw
   % and rjs are NOT the support length of the impuls response of the dual
   % wavelet filtering. Indeed, that impuls response would correspond to a
   % column of the forward wavelet transform. We are after the rows of that
   % wavelet transform matrix. We want to know which entries into the original
   % data set have influence on a particular wavelet coefficient.
   % Probably the easiest way to think about this, is by reasoning on the
   % non-decimated transform. Since wavelet coefficients of that non-decimated
   % transform are interlaced versions (cycle spinning) of shifted fast wavelet
   % transforms, the coefficients in that redundant transform depend on the
   % same number of input values as the fast wavelet coefficients at the same
   % scale. In order to compensate for the downsampling, the filters are
   % upsampled in each step, which explains why their length increases as
   % 2^(J-1-j)*ngt-2^(J-1-j)+1.
   firsts = first(mod((1:2:2^(j+1))+sht-1,2^(j+1))+1);
   if rjs > N,
      lasts = mod(firsts-2,N)+1;
   else
      lasts = last(mod((1:2:2^(j+1))+nht+sht-2,2^(j+1))+1);
   end
   firstw = first(mod((1:2:2^(j+1))+sgt-1,2^(j+1))+1);
   if rjw > N,
      lastw = mod(firstw-2,N)+1;
   else
      lastw = last(mod((1:2:2^(j+1))+ngt+sgt-2,2^(j+1))+1);
   end

   first(1:2^j) = firsts;
   first(2^j+1:2^(j+1)) = firstw;
   last(1:2^j) = lasts;
   last(2^j+1:2^(j+1)) = lastw;
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
