
function [wt, thr] = minmsehard(v,w)

% minmsehard -- apply hard threshold with minimum MSE threshold
%                    MSE = Loss = mean squared error
%  Usage
%    [wt, thr] = minmsehard(v,w);
%  Inputs
%    v      noise-free input vector
%    w      observed input vector
%  Outputs
%    wt     thresholded vector
%    thr    minimum MSE threshold
%  Description
%    Applies minimum MSE hard threshold to all elements of w
%  Note
%
%  See also
%    help MSEhard
%    help minmsethreshest

[wt, thr] = minmsethreshest(v,w,0);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 

