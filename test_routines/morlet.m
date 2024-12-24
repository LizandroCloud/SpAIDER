% Set effective support and grid parameters. 
lb = -4; ub = 4; n = 1000; 
% Compute and plot Morlet wavelet. 
[psi,x] = morlet(lb,ub,n); 
plot(x,psi), title('Morlet wavelet')