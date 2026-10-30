
function [y, thr] = minsure(x,stdev)

% minsure -- soft thresholding with min sure threshold in 1D
% Obsolete routine. Use [y, thr] = minsuresoft(x,stdev) instead

if nargin < 2,
   stdev = 1;
end

[y, thr] = minsuresoft(x,stdev);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
