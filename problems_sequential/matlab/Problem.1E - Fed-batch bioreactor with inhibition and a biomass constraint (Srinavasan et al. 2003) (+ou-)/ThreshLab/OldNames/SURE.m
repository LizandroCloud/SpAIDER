
function sure = SURE(x,t,stdev)

% SURE -- apply soft threshold + compute Stein's unbiased risk estimator (SURE)
% Obsolete routine. Use sure = SUREsoft(x,t,stdev); instead

if nargin == 2,
   stdev = 1;
end

sure = SURE1soft(x,t,stdev);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
