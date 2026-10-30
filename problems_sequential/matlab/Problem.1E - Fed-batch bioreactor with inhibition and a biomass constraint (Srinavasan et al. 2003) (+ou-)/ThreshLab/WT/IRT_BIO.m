
function f = IRT_BIO(w,h,g,sh,sg,alfa)

% usage: f = IRT_BIO(w,h,g,sh,sg,alfa);
% or:    f = IRT_BIO(w,h,g,sh,sg);
% 
% Inverse wavelet transform
% Non-decimated wavelet transform, redundant wavelet transform, stationary
% wavelet transform
% input is matrix with dimensions (M+1) x N
% first M rows contain wavelet coefficients
% row (M+1) contains scaling coefficients
% output is row vector with length N
% non-decimated data are redundant => inverse transform is not unique.
% Paramater alfa controls which reconstruction is taken. Use alfa = 0.5
% for first application and study the code if you want to explore
% (alfa = 0.5 is default value)

if nargin == 5,
   alfa = 0.5;
end

M = size(w,1)-1;

s = w(M+1,:,:);
N = size(s,2);
ph = length(h);
pg = length(g);

%%% SPLITTING THE FILTERS %%%
heven = h(1:2:ph);
hodd  = h(2:2:ph);
geven = g(1:2:pg);
godd  = g(2:2:pg);
   
p = div(max(ph,pg)+1,2)*2;
% p must be even for two reasons: p determines how much an input is extended.
% If this extension is not even, vectors of nonsimilar lengths are being added
% and also we are adding a shift of p/2

skip = 2^(M-1);
for j = M:-1:1
  s = shift(w(M+1,:),-skip*sh);
  wj = shift(w(j,:),-skip*sg);
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
      sk = s(k:2*skip:Njp);
      wjk = wj(k:2*skip:Njp);
      x1(k:2*skip:Njp) = filter(heven,1,sk) + filter(geven,1,wjk);
      x1(k+skip:2*skip:Njp) = filter(hodd,1,sk) + filter(godd,1,wjk);
      sk = s(skip+k:2*skip:Njp);
      wjk = wj(skip+k:2*skip:Njp);
      x2(k:2*skip:Njp) = filter(heven,1,sk) + filter(geven,1,wjk);
      x2(k+skip:2*skip:Njp) = filter(hodd,1,sk) + filter(godd,1,wjk);
  end
  x1(:,1:skip*p) = x1(:,Nj+1:Nj+skip*p);
  x2(:,1:skip*p) = x2(:,Nj+1:Nj+skip*p);
  x1 = x1(1:N);
  x2 = shift(x2(1:N),-skip);
  x = alfa*x1 + (1-alfa)*x2;
  w(M+1,:) = x(1:N);
  skip = skip/2;
end
f = w(M+1,:);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
