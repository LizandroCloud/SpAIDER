wname = 'db1';
[phi,psi,xval] = wavefun(wname,7);
subplot(211); stairs(xval,psi); title('Wavelet'); 

% Compute and plot wavelet integrals approximations 
% on a dyadic grid. 
[integ,xval] = intwave(wname,21); 
subplot(212); stairs(xval,integ); 
title(['Wavelet integrals over [-Inf x] ' ... 
       'for each value of xval']);