
function f = condprob1(X,W,sig)

% approximated conditional probabilty (without normalizing): mu is considered
% to be infinity.
% see condprob.m


f1 = 1-cumgauss(W+sig,sig)+cumgauss(W-sig,sig);
% f1 has to be normalized

f0 = (cumgauss(W+sig,sig)-cumgauss(W-sig,sig))/(2*sig);

f  = X.*f1 + (1-X).*f0;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
