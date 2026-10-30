
function pl = PLising(X,tau);

% pl = PLising(X,tau);
% pseudo-likelihood function for Ising model with parameter tau
% X is binary image with pixel values 0 and 1

Uc0 = isingpotupdate(2*X-1,2*X-1,tau);
Uc1 = isingpotupdate(2*X-1,1-2*X,tau);
localprobs = exp(-Uc0)./( exp(-Uc0) + exp(-Uc1) );
pl = sum(sum(log(localprobs)));

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
