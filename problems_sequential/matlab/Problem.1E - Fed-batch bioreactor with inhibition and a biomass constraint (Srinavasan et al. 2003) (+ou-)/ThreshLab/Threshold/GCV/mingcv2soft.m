
function [y, thr] = mingcv2soft(x)

% mingcv2soft -- minimizes gcv in 2D
%  Usage
%    [y, thr] = mingcv2soft(x)
%  Inputs
%    x      input coefficients
%  Outputs
%    y      the shrunk output coefficients
%  Description
%    x is shrunk with threshold thr and placed in y
%    inverse wavelet transform of y gives result
%  Note
%    The optimal thr is found by Fibonacci search for gcv function
%    This minimization procedure is only a 'rough' one

y = x;

figure(1)
%initialisation Fibonacci - numbers
F0 = 1;
F1 = 1;
F2 = 2;
Niter = 1;

m = 50;
b = max(max(abs(x)));
a =  b / m;
eps = 0.0001;
while (b - a) / F2 > eps
  F0 = F1;
  F1 = F2;
  F2 = F0 + F1;
  Niter = Niter + 1;
end

v = a + (F1/F2) * (b - a);
fv = GCV2(x,v);
u = a + (F0/F2) * (b - a);
fu = GCV2(x,u);

for k = 1:Niter - 1
 if (fu > fv)
    a = u;
    u = v;
    fu = fv;
    v = b-u+a;
    fv = GCV2(x,v);
 else
    b = v;
    v = u;
    fv = fu;
    u = a+b-v;
    fu = GCV2(x,u);
 end
end

thr = b;


y = ST(x,thr);


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
