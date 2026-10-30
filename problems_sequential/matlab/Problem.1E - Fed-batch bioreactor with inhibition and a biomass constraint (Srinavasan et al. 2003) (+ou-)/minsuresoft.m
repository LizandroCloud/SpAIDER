
function [y, thr] = minsuresoft(x,stdev)

% minsuresoft -- soft thresholding with min sure threshold in 1D
%  Usage
%    [y, thr] = minsuresoft(x,stdev);
%    y = minsuresoft(x,stdev);
%    [y, thr] = minsuresoft(x);
%    y = minsuresoft(x);
%  Inputs
%    x      input coefficients
%    stdev  standard deviation of the noise (or estimation)
%           default value is 1
%  Outputs
%    y      the shrunk output coefficients
%    thr    threshold minimizing SUREsoft(thr)
%  Description
%    x is shrunk with threshold thr and placed in y
%    Inverse wavelet transform of y gives result
%  Note
%    The optimal thr is found by Fibonacci search for SURE function
%  See also
%    help mingcvsoft
%    help SURE

if nargin < 2,
   stdev = 1;
end

y = x;

%initialisation Fibonacci numbers
F0 = 1;
F1 = 1;
F2 = 2;
Niter = 1;

% computation of number of iteration steps in minimisation
m = 50;
b = max(abs(x)); 
a =  b / m; %% 
eps = 0.0001;
while (b - a) / F2 > eps
  F0 = F1;
  F1 = F2;
  F2 = F0 + F1;
  Niter = Niter + 1;
end

firstmove = 0;
while ~firstmove

     v = a + (F1/F2) * (b - a);
     fv = SUREsoft(x,v,stdev); 
     u = a + (F0/F2) * (b - a);
     fu = SUREsoft(x,u,stdev);
     
     for k = 1:Niter - 1
      if (fu > fv)
         a = u;
         u = v;
         fu = fv;
         v = b-u+a;
         fv = SUREsoft(x,v,stdev);
     firstmove = 1;
      else
         b = v;
         v = u;
         fv = fu;
         u = a+b-v;
         fu = SUREsoft(x,u,stdev);
      end
     end
     
     if (b < 10^(-10))
       firstmove = 1;
     end % if
     
     a = b/m;
     
end %while ~firstmove


thr = b;


y = ST(x,thr);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
