
clear all
N = 2048;
f = makesignal((1:N),3);

signal = f;
figure(1)
subplot(1,1,1)
plot(signal)
axis([0,N,-2,12])
title('Uncorrupted signal')
axis([0,N,-2,12])
