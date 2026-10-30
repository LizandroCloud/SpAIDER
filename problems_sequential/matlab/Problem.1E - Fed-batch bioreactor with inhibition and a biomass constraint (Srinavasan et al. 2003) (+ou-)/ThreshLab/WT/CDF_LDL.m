function [h, ht] = CDF_LDL(n,nt)

% usage: [h, ht] = CDF_LDL(n,nt);
% this function is called by makeWTfilter.m 
% Best is to use makeWTfilter.m 
% creates Cohen-Daubechies-Feauveau scaling filters `with less dissimilar
% lengths'
% !!!!!!!!!! at this moment only the 4-4 version is available
% n is number of primal vanishing moments.
% nt is number of dual vanishing moments.
% n+nt most be even

nn = n + nt;

if rem(nn,2) 
  error('n + nt must be even');
end
nn = nn/2;

ht = [0.026749 -0.016864 -0.078223 0.266864 0.602949 0.266864 -0.078223 ...
 -0.016864 0.026749];
ht = ht*sqrt(2);
h = [-0.045636 -0.028772 0.295636 0.557543 0.295636 -0.028772 -0.045636];
h = h*sqrt(2);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
