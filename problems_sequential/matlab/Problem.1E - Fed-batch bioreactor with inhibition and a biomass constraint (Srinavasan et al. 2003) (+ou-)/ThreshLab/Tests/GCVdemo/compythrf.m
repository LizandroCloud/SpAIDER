
subplot(1,1,1)
plot(signal)
title('Comparison thresholding output with uncorrupted signal.')
hold on
plot(f,'r')
hold off
axis([0,2048,-2,12])
