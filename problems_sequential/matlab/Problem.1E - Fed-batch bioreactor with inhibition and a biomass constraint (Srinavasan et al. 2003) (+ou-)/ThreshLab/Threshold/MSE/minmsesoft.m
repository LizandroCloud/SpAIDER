
function [wt, thr] = minmsesoft(v,w)

% minmsesoft -- apply soft threshold with minimum MSE threshold
%                    MSE = Loss = mean squared error
%  Usage
%    [wt, thr] = minmsesoft(v,w);
%  Inputs
%    v      noise-free input vector
%    w      observed input vector
%  Outputs
%    wt     thresholded vector
%    thr    minimum MSE threshold
%  Description
%    Applies minimum MSE soft threshold to all elements of w
%  Note
%
%  See also
%    help MSEsoft
%    help minmsethreshest

[wt, thr] = minmsethreshest(v,w,1);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 

