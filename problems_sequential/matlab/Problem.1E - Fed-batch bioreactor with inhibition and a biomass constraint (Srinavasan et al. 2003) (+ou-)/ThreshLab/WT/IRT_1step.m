
function x = IRT_1step(s,w,j,h,g,sh,sg,alfa)

% v = IRT_1step(s,w,j,h,g,sh,sg,alfa)
% 
% Inverse wavelet transform, one step
% s: scaling coefficients at scale j
% w: wavelet coefficients at scale j
% j input resolution level must be a positive non-zero integer (j = 0
% corresponds to resolution of non-transformed data)
% Non-decimated wavelet transform, redundant wavelet transform, stationary
% wavelet transform

N = size(s,2);
ph = length(h);
pg = length(g);

%%% SPLITTING THE FILTERS %%%
heven = h(1:2:ph);
hodd  = h(2:2:ph);
geven = g(1:2:pg);
godd  = g(2:2:pg);
   
p = div(max(ph,pg)+1,2)*2;
% p must be even for two reasons:
% otherwise vectors of nonsimilar lengths are being added)
% and also we are adding a shift of p/2

skip = 2^(j-1);
s = shift(s,-skip*sh);
wj = shift(w,-skip*sg);
% periodic exetension
while (p*2*skip > size(s,2))
    s = [s s(:,1:N)];
    wj = [wj wj(:,1:N)];
end
Nj = size(s,2);
extra = 2*skip-mod(Nj,2*skip);
s = [s s(:,1:extra)];
wj = [wj wj(:,1:extra)];
s = [s s(:,extra+1:extra+2*skip*p)];
wj = [wj wj(:,extra+1:extra+2*skip*p)];
x1 = s;
x2 = s;
Njp = size(wj,2);
for k = 1:skip
    sk = s(:,k:2*skip:Njp);
    wjk = wj(:,k:2*skip:Njp);
    x1(:,k:2*skip:Njp) = filter(heven,1,sk')' + filter(geven,1,wjk')';
    x1(:,k+skip:2*skip:Njp) = filter(hodd,1,sk')' + filter(godd,1,wjk')';
    sk = s(:,skip+k:2*skip:Njp);
    wjk = wj(:,skip+k:2*skip:Njp);
    x2(:,k:2*skip:Njp) = filter(heven,1,sk')' + filter(geven,1,wjk')';
    x2(:,k+skip:2*skip:Njp) = filter(hodd,1,sk')' + filter(godd,1,wjk')';
end
x1(:,1:skip*p) = x1(:,Nj+1:Nj+skip*p);
x2(:,1:skip*p) = x2(:,Nj+1:Nj+skip*p);
x1 = x1(:,1:N);
x2 = shift(x2(:,1:N),-skip);
x = alfa*x1 + (1-alfa)*x2;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
