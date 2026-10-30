
function x = IFHFW(w,L,h,g,sh,sg);

% usage: x = IFHFW(w,L,h,g,sh,sg);
%
% Inverse Fast Haar-Fisz-Wavelet Transform, as proposed in
% Fryzlewicz, P. and Nason, G. P., A Haar-Fisz algorithm for Poisson intensity 
% estimation. Journal of Computational and Graphical Statistics.
%
% See also: FFHFW.m

xhf = IWT_BIO(w,L,h,g,sh,sg);
hh = [1 1]/sqrt(2);
gg = [1 -1]/sqrt(2);
whf = FWT_BIO(xhf,0,hh,gg,0,0);
x = IFMF(whf,0);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
