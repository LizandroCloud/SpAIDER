
figure(1)

if exist('n_gcv_points') ~= 1, n_gcv_points = 50; end
m = n_gcv_points;
thmax = max(abs(waveletcoef(2^L+1:N)));
thmax = thmax/8;

disp('Compute and estimate output error');
disp('...in practical applications, exact error computation is impossible');
disp('   since the uncorrupted signal is unknown...');
th = (1:m) * thmax / m;
wthr = waveletcoef;
uncorruptedcoef = FWT_BIO(f,L,ht,gt,sht,sgt);
gcv = GCV1(waveletcoef((2^(L)+1):N),th);
mse = MSE1(uncorruptedcoef,waveletcoef,th);
subplot(1,1,1)
plot(th,gcv,'r')
hold on
plot(th,mse)
hold off
title('looking for minimum error threshold')
text(th(fix(m*0.70)),gcv(fix(m*0.75))*1.1,'estimated error (GCV)')
text(th(fix(m*0.70)),mse(fix(m*0.75))*1.1,'exact error')

threshold = mingcvthresh(waveletcoef(2^L+1:N));
GCVthreshold = threshold
