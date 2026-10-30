
function wcoef = RWT_1step(x,j,gt,sgt)

% wcoef = RWT_1step(x,j,gt,sgt)
% wcoef = RWT_1step(x,j,ht,sht)
% 
% Non-decimated wavelet transform, redundant wavelet transform, stationary
% wavelet transform: THIS PROCEDURE PERFORMS ONE STEP (from scale j to scale
% j+1, where scale j = 1 corresponds to the finest resolution.) and only one
% subband at the time (lowpass OR highpass)
% When x is a matrix, a wavelet transform is performed on all rows of x.
% See RWT_BIO

N = size(x,2);
Nj = N;
pgt = length(gt);
sgt = -sgt - pgt + 1;
   
% filters must be flipped in forward wavelet transform
% flipping causes change in shifts sht and sgt (see above)
gt = gt(pgt:-1:1);
p = pgt + mod(pgt,2);
% p must be even for two reasons:
% otherwise vectors of nonsimilar lengths are being added)
% and also we are adding a shift of p/2

wcoef = zeros(size(x));

skip = 2^(j-1);
x = shift(x,-skip*sgt);
% periodic exetension
while (p*skip > size(x,2))
    for tel = 1:skip
       x = [x x(:,1:N)];
    end
Nj = size(x,2);
end
x = [x x(:,1:skip*p)];
wj = zeros(size(x));
Njp = size(wj,2);
for k = 1:skip
    xk = x(:,k:skip:Njp);
    wj(:,k:skip:Njp) = filter(gt,1,xk')';
end
wj(:,1:skip*p) = wj(:,Nj+1:Nj+skip*p);
wcoef = wj(:,1:N);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
