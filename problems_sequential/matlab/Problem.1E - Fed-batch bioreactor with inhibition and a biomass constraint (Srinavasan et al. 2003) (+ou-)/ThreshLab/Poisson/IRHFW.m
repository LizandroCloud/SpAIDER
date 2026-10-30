
function x = IRHFW(w,h,g,sh,sg,alfa);

% usage: x = IRHFW(w,h,g,sh,sg,alfa);
%
% Inverse redundant Haar-Fisz-Wavelet Transform, as proposed in
% Fryzlewicz, P. and Nason, G. P., A Haar-Fisz algorithm for Poisson intensity 
% estimation. Journal of Computational and Graphical Statistics.
%
% see also: FRHFW.m

if nargin == 5,
   alfa = 0.5;
end

M = size(w,1)-1;
J = log2(size(w,2));
L = J-M;
xhf = IRT_BIO(w,h,g,sh,sg,alfa);
hh = [1 1]/sqrt(2);
gg = [1 -1]/sqrt(2);
whf = FWT_BIO(xhf,0,hh,gg,0,0);
x = IFMF(whf,0);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
