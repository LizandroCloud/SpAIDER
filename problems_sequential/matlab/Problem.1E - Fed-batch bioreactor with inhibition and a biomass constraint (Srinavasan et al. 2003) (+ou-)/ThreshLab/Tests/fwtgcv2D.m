N = size(input,1);
J = round(log(N)/log(2));
L = J - M;
wcoef = FWT2_BIO(input,L,ht,gt,sht,sgt);
wout = wcoef;
clear gcvthr

for j = L:J-1
% L is coarsest level; J-1 is finest level
% take L = 5: number of coefficients at coarsest level is then
% 2^L*2^L = 1024, which is about the minimum for successful application of GCV
% (asymptotic method!)
j

  % horizontal component

  x = wcoef(2^j+1:2^(j+1),1:2^j);
  [y, thr] = mingcv2(x);
  wout(2^j+1:2^(j+1),1:2^j) = y;
  gcvthr(j,1) = thr;

  % diagonal subband

  x = wcoef(2^j+1:2^(j+1),2^j+1:2^(j+1));
  [y, thr] = mingcv2(x);
  wout(2^j+1:2^(j+1),2^j+1:2^(j+1)) = y;
  gcvthr(j,2) = thr;

  % vertical subband

  x = wcoef(1:2^j,2^j+1:2^(j+1));
  [y, thr] = mingcv2(x);
  wout(1:2^j,2^j+1:2^(j+1)) = y;
  gcvthr(j,3) = thr;

end
output = IWT2_BIO(wout,L,h,g,sh,sg);
image(output)
