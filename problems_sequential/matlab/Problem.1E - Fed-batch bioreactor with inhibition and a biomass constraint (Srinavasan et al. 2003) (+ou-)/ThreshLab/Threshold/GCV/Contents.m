% THRESHLAB Threshold/GCV
% -----------------------
% 
% Threshold selection by SURE (Stein's Unbiased Risk Estimation)
%                     and GCV (Generalized Cross Validation)
% 
% .m - files
% ----------
%    ** ST = Soft Threshold ** HT = Hard Threshold **
% Contents.m     : this file
% GCVsoft.m      : ST GCV for scalar threshold, 1D array of coefficients
% GCV1soft.m     : ST GCV in vector of thresholds (1D array of coeff.)
% GCV2soft.m     : ST GCV for scalar threshold, 2D array of coefficients
%
% GCVhard.m      : HT GCV for scalar threshold, 1D array of coefficients
% GCV1hard.m     : HT GCV in vector of thresholds (1D array of coeff.) 
%
% GCV1softhard.m : HT with ST GCV (for theoretical use only)
%
% varestGCV.m    : GCV variance estimator for data with outliers
%
% SUREsoft.m     : ST SURE for scalar threshold and 1D array of coeff.
% SURE1soft.m    : ST SURE in vector of thresholds (1D array of coeff.)
% SURE2soft.m    : ST SURE for scalar threshold, 2D array of coefficients
% 
% SUREhard.m     : ST SURE for scalar threshold and 1D array of coeff.
% SURE1hard.m    : ST SURE in vector of thresholds (1D array of coeff.)
%
% mingcvsoft.m       : minimization of GCVsoft.m
% mingcvhard.m       : minimization of GCVhard.m
% mingcv2soft.m      : minimization of GCV2soft.m
% mingcvthreshsoft.m : variant of mingcv.m (see help mingcvthresh for details)
% minsuresoft.m      : minimization of SUREsoft.m
% 
% subdirectories
% --------------
% <none>

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
