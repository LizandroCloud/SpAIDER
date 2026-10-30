
function mse = MSE(v,w,t,softORhard)

% Obsolete routine. Use mse = MSE1thresh(v,w,t,softORhard) instead

if nargin < 4, softORhard = 1; end

mse = MSE1thresh(v,w,t,softORhard);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
