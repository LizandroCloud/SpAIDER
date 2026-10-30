
clear
batse551data
x = row(x);
yp0 = row(yp0);
t = x/10^6;

N0 = length(yp0);
if ~(exist('L') == 1), L = 5; end
if length(L)>1, L = 5; end
L0 = L;
[yp L] = dyadicextension(yp0,L0);
N = length(yp);

[h, ht] = makeWTfilter('symmlet',8);
M = log2(N) - L;

yshr0 = BRCVS(yp,M,h,ht); % pilot estimator

% update pilot estimator
backgroundlevel = median(yp0);
M = ceil(log2(N0))-L0;
[h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);
w = RWT_BIO(yp,M,ht,gt,sht,sgt);
wshr = RWT_BIO(yshr0,M,ht,gt,sht,sgt);
for j = 1:M,
   j
   n = length(w(j,:));
   for i = 1:n,
      if poissoncpf(yshr0(i),backgroundlevel) > 0.95 & j > 3,
         wshr(j,i) = w(j,i);
      end
   end
end
yshr = IRT_BIO(wshr,h,g,sh,sg);


N = N0;
yp = yp(N+1:2*N);
yshr = yshr(N+1:2*N);
yshr0 = yshr0(N+1:2*N);
L = L0;


plot(t,yp,'.')
hold on
plot(t,yshr,'r')
hold off
axis([0,0.47,0,50])
figure(1)
