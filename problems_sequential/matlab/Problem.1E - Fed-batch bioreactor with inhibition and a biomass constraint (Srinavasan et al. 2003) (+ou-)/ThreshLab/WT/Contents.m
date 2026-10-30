% THRESHLAB WT
% ------------
% 
% Wavelet transform routines
% 
% .m - files
% ----------
% Contents.m     : this file
% CDFspline.m    : makes biorthogonal filters of Cohen-Daubechies-Feauveau
% CDF_LDL.m      : makes biorthogonal CDF filters with `less dissimilar lengths'
% Haar.m         : Haar filter (orthogonal)
% Daubechies.m   : Daubechies' orthogonal filters
% symmlet.m      : Daubechies' orthogonal symmlet filters
% coiflet.m      : Daubechies' orthogonal coiflet filters
% makeWTfilter.m : creates scaling functions filters
% mirrorshifts.m : creates wavelet functions filters from scaling filters
% FWT_BIO.m      : Fast Wavelet Transform; (bi)orthogonal, periodic
% IWT_BIO.m      : Inverse Wavelet Transform
% FWT2_BIO.m     : Fast 2D wavelet transform (calls FWT2_1step.m)
% FWT2_1step.m   : one step of FWT2_BIO.m (calls FWT2_halfstep.m)
% FWT2_halfstep.m: half step of FWT2_BIO.m
% IWT2_BIO.m     : Inverse 2D Wavelet Transform (calls IWT2_1step.m)
% IWT2_1step.m   : one step of IWT2_BIO.m (calls IWT2_halfstep.m)
% IWT2_halfstep.m: half step of IWT2_BIO.m
% RWT_BIO.m      : Non-Decimated Wavelet Transform 
% IRT_BIO.m      : Inverse Non-Decimated Wavelet Transform
% RWT_1step.m    : one step of RWT_BIO.m
% IRT_1step.m    : one step of IRT_BIO.m
% RWT2_BIO.m     : Non-Decimated 2D Wavelet Transform (calls RWT_1step.m)
% IRT2_BIO.m     : Inverse Non-Decimated 2D WT (calls IRT_1step.m)
% 
% plotscales.m   : plot wavelet coefficients with dashed lines at scale borders
% plotasspikes.m : plot elements of sparse vectors as spikes
% plotasblocks.m : plot as piecewise constant, rather than as polyline
%
% FWTdualsupport.m : computes support of the dual wavelet basis functions
% WTbasisfunction.m: plots mother or father function at given resolution
% 
% 
% subdirectories
% --------------
% <none>

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
