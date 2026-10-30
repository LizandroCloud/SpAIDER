
N = 2048;
x = makegrid(N,0);
f = makesignal(x,3);
B = 1;
random = 0;
ratio = 10;
y = makenoisy(f,B,random,ratio);
[h ht] = Daubechies(3);
[h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);
J = log2(N);
M = 5;
L = J - M;
w = FWT_BIO(y,L,ht,gt,sht,sgt);
v = FWT_BIO(f,L,ht,gt,sht,sgt);

wthr = w;
wthr(2^L+1:N) = mingcv(w(2^L+1:N));
ythr = IWT_BIO(wthr,L,h,g,sh,sg);

m = 100;
thrgcv = mingcvthresh(w(2^L+1:N));
thrmax = 5*thrgcv;
thr = (1:m)/m * thrmax;
gcv = GCV1(w(2^L+1:N),thr);
mse = MSE1(v(2^L+1:N),w(2^L+1:N),thr);

figure(2)
subplot(1,1,1)
plot(thr,gcv,'--')
hold on
plot(thr,mse)
hold off
legend('GCV','MSE')


figure(1)
subplot(1,1,1)
plot(y)
hold on
plot(f,'r')
hold off
axis([0,N,-2,12])


figure(3)
subplot(1,1,1)
plot(ythr)
ratio1 = snr(f,y-f)
ratio2 = snr(f,ythr-f)
axis([0,N,-2,12])
