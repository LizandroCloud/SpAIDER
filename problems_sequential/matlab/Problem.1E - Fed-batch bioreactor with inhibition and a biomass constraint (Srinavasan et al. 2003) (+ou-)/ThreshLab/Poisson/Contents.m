% THRESHLAB Poisson
% -----------------
% 
% Routines for Poisson intensity estimation
% 
% .m - files
% ----------
% Contents.m     : this file
% poisson.m      : Poisson random number generator
% FWTpoisnorm.m  : Conditional Variance Stabilisation for Fast WT.
% RWTpoisnorm.m  : Conditional Variance Stabilisation for Redundant WT.
% EsqrtPoisson.m : computes E(sqrt(X+c)), where X ~ Poisson(lambda)
% VsqrtPoisson.m : computes V(sqrt(X)), where X ~ Poisson(lambda)
% VsqrtPoissonN0.m : computes V(sqrt(X) | X > 0), where X ~ Poisson(lambda)
% FFMF.m         : Forward Fast Multiscale Fisz Transform
% IFMF.m         : Inverse Fast Multiscale Fisz Transform
% FRMF.m         : Forward Redundant (= non-decimated) Multiscale Fisz Tr.
% IRMF.m         : Inverse Redundant (= non-decimated) Multiscale Fisz Tr.
% FFHFW.m        : Forward Fast Haar-Fisz-Wavelet Transform
% IFHFW.m        : Inverse Fast Haar-Fisz-Wavelet Transform
% FRHFW.m        : Forward Redundant (= non-decimated) Haar-Fisz-Wavelet Tr.
% IRHFW.m        : Inverse Redundant (= non-decimated) Haar-Fisz-Wavelet Tr.
% medianmass0norm: Median of r.v. with mixture density: point mass + normal
% CDFmass0norm   : Cum. Distr. Fun. of r.v. with mixture density
% BFCVS.m        : Bayesian Fast (= decimated) Conditional Var. Stabilisation
% BRCVS.m        : Bayesian Redundant Conditional Variance Stabilisation
% MmLEmixnorm    : Max. marginal likelihood est. of mixture parameter in
%                  mixed normal model
% lnLmixednormal.m: Log-likelikihood as function of mixture parameter in
%                   mixed normal model
% dlnLmixednormaldp: Derivative of lnLmixednormal(p)
% 
% subdirectories
% --------------
% PoissonTests

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
