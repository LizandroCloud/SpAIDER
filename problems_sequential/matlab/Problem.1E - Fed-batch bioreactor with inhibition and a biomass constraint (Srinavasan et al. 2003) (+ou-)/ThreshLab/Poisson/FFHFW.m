
function w = FFHFW(x,L,ht,gt,sht,sgt)

% usage: w = FFHFW(x,L,ht,gt,sht,sgt);
%
% Forward Fast Haar-Fisz-Wavelet Transform, as proposed in
% Fryzlewicz, P. and Nason, G. P., A Haar-Fisz algorithm for Poisson intensity 
% estimation. Journal of Computational and Graphical Statistics.

% First, use a fast multiscale Fisz transform for variance stabilisation and
% gaussianisation. 
% Since Multiscale Fisz is just used as a variance stabilisation, the
% decomposition depth has to be maximal, i.e., up to a single scaling function
% at the coarsest possible level
whf = FFMF(x,0);
% whf is homoscedastic and close to uncorrelated. Use a normalised inverse Haar
% transform + normalised wavelet decomposition to maintain this
% quasi-stationarity
hh = [1 1]/sqrt(2);
gg = [1 -1]/sqrt(2);
xhf = IWT_BIO(whf,0,hh,gg,0,0);
w = FWT_BIO(xhf,L,ht,gt,sht,sgt);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
