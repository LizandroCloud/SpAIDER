% Set number of iterations and wavelet name. 
iter = 10;
wav = 'db7';
name = 'Daubechies';
% Compute approximations of the wavelet function using the
% cascade algorithm. 
for i = 1:iter 
    [phi,psi,xval] = wavefun(wav,i); 
%     plot(xval,psi); 
%     hold on 
%     plot(xval,phi); 
end
fig=figure(1);
axes1 = axes('Parent',fig);
FigHandle = fig;
set(FigHandle, 'Position', [100, 100, 500, 300]);
plot1 = plot(xval,psi,'k','Parent',axes1);
set(plot1(1),'DisplayName','Wavelet mãe');
hold on;
plot2 = plot(xval,phi,'--k','Parent',axes1);
set(plot2(1),'DisplayName','Wavelet escala');

title(['Wavelets  ',name]); 
legend(axes1,'show');
hold on;
figure(10)
[wp,x] = wpfun(wav,7);

% 
% % Plot wavelet function. 
% [phi,psi,xval] = wavefun(wname,7);
% subplot(211); plot(xval,psi); title('Wavelet'); 
% 
% % Compute and plot wavelet integrals approximations 
% % on a dyadic grid. 
% [integ,xval] = intwave(wname,7); 
% subplot(212); plot(xval,integ); 
% title(['Wavelet integrals over [-Inf x] ' ... 
%        'for each value of xval']);