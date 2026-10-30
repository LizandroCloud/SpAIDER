
function [wt thr] = oracleselection(v,w,stdev)

% oracleselection -- ThreshLab/Threshold/ - oracle variable selection from 
%                    observations in model w = v + error 
%  Usage
%    [wt, thr] = oracleselection(v,w);
%    [wt, thr] = oracleselection(v,w,stdev);
%  Inputs
%    v      error-free input vector
%    w      observed input vector
%    stdev  standard deviation of errors w-v
%  Outputs
%    wt     vector with size equal to input
%    thr    minimum MSE threshold
%  Description
%    selects components from w if error-free value is smaller than stdev.
%    The information whether |v|< stdev, is of course "secret", i.e., hidden in
%    practical situations. We need an "oracle" to know this.
%    stdev acts as a "threshold" on v (not on w)
%  Note
%    The routine does not assume full knowledge of v
%  See also
%    help minmsehard

if nargin < 3,
   stdev = sqrt(var2(w-v));
end
wt = w.*(abs(v)>stdev);
thr = stdev;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
