
function [wt, thr] = minmsethreshest(v,w,softORhard,varargin)

% minmsethreshest -- apply soft/hard threshold with minimum MSE threshold
%                    MSE = Loss = mean squared error
%  Usage
%    [wt, thr] = minmsethreshest(v,w,softORhard,'maxthr',maxthr,'minthr',minthr);
%  Inputs
%    v      noise-free input vector
%    w      observed input vector
%    softORhard boolean: 1 for soft; 0 for hard
%    maxthr maximum possible threshold (optional)
%    minthr minimum possible threshold (optional)
%  Outputs
%    wt     thresholded vector
%    thr    minimum MSE threshold
%  Description
%    Applies minimum MSE soft- or hard threshold to all elements of w 
%  Note
%    For making a plot of mse for a given vector of thresholds
%    use MSE1thresh (for one dimension)
%
%  See also
%    MSEthresh

if nargin<3, softORhard = 1; end

N = size(w,1)*size(w,2);
wt = w;

% identification of maximum and minimum threshold
% default values:
maxthr = max(abs(reshape(w,N,1)));
minthr = 0;
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
n_iter = 1;
a = minthr;
b = maxthr;

% Computation of needed number of iteration steps
eps = 0.0001;
while (b - a) / F2 > eps
  F0 = F1;
  F1 = F2;
  F2 = F0 + F1;
  n_iter = n_iter + 1;
end

% minimization routine

r = a + (F1/F2) * (b - a);
fr = MSEthresh(v,w,r,softORhard); 
l = a + (F0/F2) * (b - a);
fl = MSEthresh(v,w,l,softORhard);
     
for k = 1:n_iter - 1
   if (fl > fr)
      a = l;
      l = r;
      fl = fr;
      r = b-l+a;
      fr = MSEthresh(v,w,r,softORhard);
   else
      b = r;
      r = l;
      fr = fl;
      l = a+b-r;
      fl = MSEthresh(v,w,l,softORhard);
   end
end
     
thr = b;

if softORhard==1,
   wt = ST(w,thr);
else
   wt = HT(w,thr);
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
