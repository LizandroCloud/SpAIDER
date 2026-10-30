
function v = IWT2_BIOhalfstep(w,h,g,sh,sg);

[m n] = size(w);
J = log2(n);
J = round(J);
v = IWT_BIO(w,J-1,h,g,sh,sg);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
