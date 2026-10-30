
wthr = waveletcoef;
N = length(waveletcoef);
wthr((2^(L)+1):N) = ST(waveletcoef((2^(L)+1):N),threshold);
figure(2)
subplot(1,1,1)
plotscales(wthr(2:length(wthr)),7,wm,wM);
title('Thresholded coefficients (soft-thresholding)')
signal = IWT_BIO(wthr,L,h,g,sh,sg);
figure(1)
subplot(1,1,1)
plot(signal)
axis([0,N,-2,12])
title('The output')
figure(1)
