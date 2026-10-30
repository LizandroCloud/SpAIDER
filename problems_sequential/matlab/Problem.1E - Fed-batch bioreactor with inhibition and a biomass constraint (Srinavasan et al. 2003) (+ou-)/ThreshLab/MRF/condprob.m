
function f = condprob(X,W,sig,mu)

% f = condprob(X,W,sig,mu);
% see PhD page 144
% original definition

f1 = (cumgauss(W+mu,sig)-cumgauss(W+sig,sig)+cumgauss(W-sig,sig)-cumgauss(W-mu,sig))/(2*(mu-sig));

f0 = (cumgauss(W+sig,sig)-cumgauss(W-sig,sig))/(2*sig);

f  = X.*f1 + (1-X).*f0;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
