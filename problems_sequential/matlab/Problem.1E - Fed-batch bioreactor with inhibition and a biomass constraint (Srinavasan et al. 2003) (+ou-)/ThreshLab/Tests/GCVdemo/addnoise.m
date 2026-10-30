
sig = 0.5;
noise = randn(size(f))*sig;
signal = f+noise;
figure(1)
subplot(1,1,1)
plot(signal)
hold on
plot(f,'r')
hold off
title('Noisy signal')
axis([0,N,-2,12]);
