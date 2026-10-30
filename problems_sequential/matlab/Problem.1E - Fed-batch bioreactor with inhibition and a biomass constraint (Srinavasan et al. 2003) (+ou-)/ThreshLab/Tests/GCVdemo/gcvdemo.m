
figure(1)
disp('make demo signal...')
makedemosignal
disp('done.')
pause(-1)
disp('computing wavelet transform uncorrupted signal (Haar)...')
WTdemo
disp('done.')
pause(-1)
disp('Add noise...')
addnoise
disp('done.')
pause(-1)
disp('wavelet transform of noisy signal...')
WTdemo
disp('done.')
ans = 1;
input('choose threshold value :');
threshold = abs(ans);
figure(2)
hold on
plot(threshold*ones(size(waveletcoef)),'r--');
plot(-threshold*ones(size(waveletcoef)),'r--');
hold off
disp('Press <Enter> or <Return> to apply this threshold.');
pause(-1)
disp('Thresholding and inverse wavelet transform...');
reducenoise
disp('ress <Enter> or <Return> for comparison with uncorrupted signal.');
pause(-1)
figure(1)
compythrf
disp('Press <Enter> or <Return> to find best (=minimum GCV) threshold');
pause(-1)
figure(1)
findbest
figure(2)
plotscales(waveletcoef(2:length(waveletcoef)),7,wm,wM);
title('wavelet coefficients')
hold on
plot(threshold*ones(size(waveletcoef)),'r--');
plot(-threshold*ones(size(waveletcoef)),'r--');
hold off
disp('Press <Enter> or <Return> to apply this threshold.');
pause(-1)
disp('Thresholding and inverse wavelet transform...');
reducenoise
disp('ress <Enter> or <Return> for comparison with uncorrupted signal.');
pause(-1)
figure(1)
compythrf
disp('THE END.')
figure(1)

