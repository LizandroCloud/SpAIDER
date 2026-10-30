
function [thr, y] = mingcvthresh(x)

% mingcv -- soft thresholding with min gcv threshold in 1D
% Obsolete routine. Use [y, thr] = mingcvthreshsoft(x) instead

[y thr] = mingcv(x);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
