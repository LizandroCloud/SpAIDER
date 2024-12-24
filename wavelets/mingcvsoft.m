
function [y, thr] = mingcvsoft(x,varargin)

% mingcvsoft -- soft thresholding with min gcv threshold in 1D
%  Usage
%    [y, thr] = mingcvsoft(x,'maxthr',maxthr,'minthr',minthr);
%    y = mingcvsoft(x);
%  Inputs
%    x      input coefficients
%    maxthr maximum possible threshold (optional)
%    minthr minimum possible threshold (optional)
%  Outputs
%    y      the shrunk output coefficients
%    thr    threshold minimizing GCVsoft(thr)
%  Description
%    x is shrunk with threshold thr and placed in y
%    Inverse wavelet transform of y gives result
%  Note
%    The optimal thr is found by Fibonacci search for gcv function
%  See also
%    help mingcvthreshsoft
%    help GCVsoft

y = x;

% identification of maximum and minimum threshold
% default values:
m = 50;
maxthr = max(abs(x));
minthr = maxthr / m; %% GCVsoft(0) is undefined
nvarargin = length(varargin);
if mod(nvarargin,2), 
   warning('inconsistent specification of variable input arguments')
end
for k = 1:div(nvarargin,2),
   varname = varargin{2*k-1};
   switch varname(1:6),
   case {'maxthr','thrmax'},
      maxthr = varargin{2*k};
      minthr = maxthr/m;
   case {'minthr','thrmin'},
      minthr = varargin{2*k};
   end
end

%initialisation Fibonacci numbers
F0 = 1;
F1 = 1;
F2 = 2;
Niter = 1;
a = minthr;
b = maxthr;

% computation of number of iteration steps in minimisation
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
     fv = GCVsoft(x,v); 
     u = a + (F0/F2) * (b - a);
     fu = GCVsoft(x,u);
     
     for k = 1:Niter - 1
      if (fu > fv)
         a = u;
         u = v;
         fu = fv;
         v = b-u+a;
         fv = GCVsoft(x,v);
     firstmove = 1;
      else
         b = v;
         v = u;
         fv = fu;
         u = a+b-v;
         fu = GCVsoft(x,u);
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
