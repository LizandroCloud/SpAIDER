
clear
makehitsplot
yp0 = row(hits);
% text(12,1400,'1')
% text(170,6500,'2')

N0 = length(yp0);
if ~(exist('L') == 1), L = 5; end
if length(L)>1, L = 5; end
L0 = L;
[yp L] = dyadicextension(yp0,L0);
N = length(yp);

[h, ht] = makeWTfilter('CDFspline',2,2);

yshr = BFCVS(yp,L,h,ht);

N = N0;
yp = yp(N+1:2*N);
yshr = yshr(N+1:2*N);
L = L0;


figure(2)
plotoption = '.';
makehitsplot
hold on
plot(yshr)
hold off
figure(1)
