function [h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);

% mirrorshifts - ThreshLab/WT/ -- finds correct detail filters for wavelet
% transforms, given scaling filters + correct filter shifts
%   Usage
%     [h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);
%     [h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h);
%   Inputs
%     h          (primal) scaling filter
%     ht         dual scaling filter
%   Outputs
%     h          (primal) scaling filter
%     g          (primal) detail filter
%     ht         dual scaling filter
%     gt         dual detail filter
%     sh         shift of primal scaling filter
%     sg         shift of primal detail filter
%     sht        shift of dual scaling filter
%     sgt        shift of dual detail filter
%   Description
%     creates all filters in forward and inverse quadrature mirror filter (QMF)
%     for a wavelet transform. Given the coefficients of the dual and primal
%     (i.e., forward and inverse/analysis and synthesis/decomposition and
%     reconstruction) two-scale equations, all other filters in the QMF and
%     their correct displacement (shifts) follow from the expressions
%     implemented in this routine.
%   Note
%     In orthogonal wavelet transforms, dual filter = primal filter, 
%     i.e., inverse transform matrix = transpose of forward matrix.
%
%     The user may want to incoroprate this procedure into  his own wavelet
%     transform routine: My_FWT = mirrorshifts + FWT_BIO
%
%   See also
%     help FWT_BIO          for the actual wavelet decomposition
%     help IWT_BIO          for the wavelet reconstruction
%     help makeWTfilter     for the design of the 2-scale equation coefficients
%                           h and ht

if (nargin == 1)
   ht = h;
end
sh  = -floor((length(h)-1)/2) ;
sht = -floor((length(ht)-1)/2) ;
clear g
p = length(h);
for i = 1:p
  gt(p-i+1) = h(i)*(-1)^(i+sh);
end
pt = length(ht);
for i = 1:pt
  g(pt-i+1) = ht(i)*(-1)^(i+sht);
end
sgt = 2 - sh - length(h);
sg = 2 - sht - length(ht);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
