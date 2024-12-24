
function wt = ST(w,t);

% usage: wt = ST(w,t);
% soft-threshold function

wt = (w > t) .* (w - t) + (w < -t) .* (w + t);


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
