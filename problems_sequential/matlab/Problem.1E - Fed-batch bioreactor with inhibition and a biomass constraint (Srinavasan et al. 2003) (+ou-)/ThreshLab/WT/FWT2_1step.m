
function ww = FWT2_onestep(x,ht,gt,sht,sgt);

% 1 step of Biorthogonal 2D, tensor-product WT

w1 = FWT2_halfstep(x,ht,gt,sht,sgt);
w2 = FWT2_halfstep(w1',ht,gt,sht,sgt);
ww = w2';

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
