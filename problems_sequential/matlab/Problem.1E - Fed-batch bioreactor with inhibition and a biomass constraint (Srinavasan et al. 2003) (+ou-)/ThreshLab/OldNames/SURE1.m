function sure = SURE1(x,t,stdev)

% SURE1 -- apply soft threshold + compute Stein's unbiased risk estimator
%                                         (SURE)
% Obsolete routine. Use sure = SURE1soft(x,t,stdev); instead

sure = SURE1soft(x,t,stdev);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
