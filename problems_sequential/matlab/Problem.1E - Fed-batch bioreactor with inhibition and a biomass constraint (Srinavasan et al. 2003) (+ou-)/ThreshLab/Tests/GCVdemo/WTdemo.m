
h = Haar;
[h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h);
L = 0;
waveletcoef = FWT_BIO(signal,L,ht,gt,sht,sgt);
figure(2)
subplot(1,1,1)
% plot(waveletcoef(2:length(waveletcoef)));
% wM = max(waveletcoef(2:length(waveletcoef)));
% wm = min(waveletcoef(2:length(waveletcoef)));
wM = 30;
wm = -30;
plotscales(waveletcoef(2:length(waveletcoef)),7,wm,wM);
title('wavelet coefficients')
