
function vv = IWT2_1step(w,h,g,sh,sg);

% 1 stap van Biorthogonale 2D, tensor-product WT

 v = IWT2_halfstep(w',h,g,sh,sg)';
vv = IWT2_halfstep(v,h,g,sh,sg);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
