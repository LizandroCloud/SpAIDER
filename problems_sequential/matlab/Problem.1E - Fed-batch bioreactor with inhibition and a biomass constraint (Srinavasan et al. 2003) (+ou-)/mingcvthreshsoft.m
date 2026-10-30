
function [thr, y] = mingcvthreshsoft(x)

% mingcvthreshsoft -- soft thresholding with min gcv threshold in 1D
%  Usage
%    [thr wthr] = mingcvthreshsoft(w);
%    thr = mingcvthreshsoft(w);
%  Inputs
%    w      input coefficients
%  Outputs
%    thr    threshold minimizing GCVsoft(thr)
%    wthr   the shrunk output coefficients
%  Description
%    minimization of GCV + soft thresholing with mingcvsoft threshold
%    inverse wavelet transform of wthr gives result
%  Note
%     identical to mingcvsoft.m, except for the order in which the output
%     appears
%  See also
%    help mingcvsoft
%    help GCVsoft

[y thr] = mingcv(x);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
