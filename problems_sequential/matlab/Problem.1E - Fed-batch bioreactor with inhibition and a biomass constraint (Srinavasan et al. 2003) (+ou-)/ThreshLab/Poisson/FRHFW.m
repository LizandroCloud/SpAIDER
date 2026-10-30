
function w = FRHFW(x,M,ht,gt,sht,sgt)

% usage: w = FRHFW(x,M,ht,gt,sht,sgt);
%
% Forward redundant Haar-Fisz-Wavelet Transform, as proposed in
% Fryzlewicz, P. and Nason, G. P., A Haar-Fisz algorithm for Poisson intensity 
% estimation. Journal of Computational and Graphical Statistics.

% First, use a _fast_ multiscale Fisz transform for variance stabilisation and
% gaussianisation. We cannot use the redundant version here: if we apply the
% inverse redundant Haar transform on redundant Fisz coefficients, we introduce
% an irriversible deformation, since the redundant Fisz coefficients are not
% the exact redudant Haar coefficients of any signal.
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
w = RWT_BIO(xhf,M,ht,gt,sht,sgt);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
