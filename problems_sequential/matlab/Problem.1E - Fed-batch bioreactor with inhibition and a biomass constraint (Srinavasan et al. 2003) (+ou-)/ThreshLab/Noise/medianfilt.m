
function y = medianfilt(x,taps)

% usage: y = medianfilt(x,taps);
% taps is integer indicating number of taps in filter

Y = row(x);
N = length(x);
M = taps;
Y = repmat(x,M,1);
m = floor(M/2);
for r=2:M
   Y(r,r:N) = x(1:N-r+1);
end
y = shift(median(Y),m);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
