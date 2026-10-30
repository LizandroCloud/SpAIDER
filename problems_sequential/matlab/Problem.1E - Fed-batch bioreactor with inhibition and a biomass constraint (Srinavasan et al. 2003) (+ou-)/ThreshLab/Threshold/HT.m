
function wt = HT(w,t);

% usage: wt = HT(w,t);
% hard-threshold function
% wt = (abs(w) > t) .* w;

wt = (abs(w) > t) .* w;


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
