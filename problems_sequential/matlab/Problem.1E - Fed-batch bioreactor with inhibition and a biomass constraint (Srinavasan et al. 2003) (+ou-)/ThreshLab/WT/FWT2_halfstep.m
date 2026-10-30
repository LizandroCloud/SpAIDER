
function w = FWT2_halfstep(x,ht,gt,sht,sgt);

% w = FWT2_halfstep(x,ht,gt,sht,sgt);

[m n] = size(x);
J = log2(n);
J = round(J);
w = FWT_BIO(x,J-1,ht,gt,sht,sgt);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
