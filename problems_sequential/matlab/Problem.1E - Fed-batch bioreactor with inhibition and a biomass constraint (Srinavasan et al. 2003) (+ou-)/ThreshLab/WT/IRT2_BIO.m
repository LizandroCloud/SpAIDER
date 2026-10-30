function f = IRT2_BIO(w,h,g,sh,sg,alfa);

% IRT2_BIO - ThreshLab/WT/ -- Inverse Redundant 2D Biorthogonal WT
%  Usage
%    f = IRT2_BIO(w,h,g,sh,sg,alfa);
%  Inputs
%    w         wavelet coefficients
%    h         Low pass sythesis filter
%    g         High pass sythesis filter
%    sh        Low pass sythesis filter shift (see below)
%    sg        High pass sythesis filter shift (see below)
%    alfa      optional design parameter (default 0.5)
%  Outputs
%    f         2D data (image/2D surface)
%  Description
%    Computes Inverse Redundant 1D Biorthogonal WT
%    The redundant wavelet transform is known under many names including
%       translation-invariant wavelet transform
%       stationary wavelet transform
%       a-trous algorithm (although this term is ambivalent)
%       cycle spinning
%       non-decimated wavelet transform
%  Note
%    sh and sg are shifts taking into account that nonzero segments of
%    filters do not necessarily start at index 1; see makeWTfilter.m and
%    mirrorshifts.m for creation of correct shifts
%
%    see RWT2_BIO.m for forward transform
%
%  See also:
%    help FWT_BIO
%    help IRT_BIO (inverse redundant transform)
%    help mirrorshifts
%    help makeWTfilter


if nargin == 5,
   alfa = 0.5;
end

M = size(w,1)-1;
n = size(w,3);
m = size(w,4);

s = zeros(n,m);
s(:,:) = w(M+1,1,:,:);
w1 = zeros(n,m);
w2 = zeros(n,m);
w3 = zeros(n,m);
for j = M:-1:1
   w1(:,:) = w(j,1,:,:);
   w2(:,:) = w(j,2,:,:);
   w3(:,:) = w(j,3,:,:);
   whor = IRT_1step(w2',w3',j,h,g,sh,sg,alfa)';
   shor = IRT_1step(s',w1',j,h,g,sh,sg,alfa)';
   s = IRT_1step(shor,whor,j,h,g,sh,sg,alfa);
end
f = s;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
