
function f = condprob2(X,W,sig,thr,p)

% f = condprob2(X,W,sig,thr,p);
% see PhD page 144
% p is probability of P( X=1 )

c = ((1-p)/p)*condprob1(0,thr,sig)/condprob1(1,thr,sig);

f1 = c*(1-cumgauss(W+sig,sig)+cumgauss(W-sig,sig));

f0 = (cumgauss(W+sig,sig)-cumgauss(W-sig,sig))/(2*sig);

f  = X.*f1 + (1-X).*f0;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
