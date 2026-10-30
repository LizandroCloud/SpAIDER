% THRESHLAB Threshold
% -----------------
% 
% Routines for noise reduction by wavelet coefficient thresholding
% 
% .m - files
% ----------
% Contents.m     : this file
% HT.m           : Hard threshold operation
% ST.m           : Soft threshold operation
% SoftreeThresh.m: Tree structured coefficient shrinking, starting from ST.
% oracle.m       : optimal coefficient selection
% 
% subdirectories
% --------------
% GCV/  : Threshold selection by SURE (Stein's Unbiased Risk Estimation)
%                             and GCV (Generalized Cross Validation)
% MSE/  : Exact mean square error plots/minimization
% Risk/ : Risk = expected MSE (theoretical issues)
% FDR/  : False Discovery Rate thresholds (Benjamini & Hochberg)
% EBayesThresh/ : Empirical Bayesian Thresholds (Johnstone & Silverman)
% IterativeThresholding/ : methods for iterative thresholding 


% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
