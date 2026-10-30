
function [Esqrt, errest, n_terms] = EsqrtPoisson(lambda,c);

% usage: Esqrt = EsqrtPoisson(lambda,c);
% computes E(sqrt(X+c)), where X ~ Poisson(lambda)
% default for c is 0

if nargin == 1,
   c = 0;
end

maxnterms = 500;
Esqrt = zeros(size(lambda));
eps = 1.e-6;
term = ones(size(lambda));
k = 0;
% computation by adding P(X=k)*sqrt(k+c) for those values in vector lambda
% where the last term was still significant. The addition is not executed for
% the other values in lambda, since this might cause NaN's. In order to prevent
% the iterarion from stopping before even the modus of the probability mass
% function has been reached, we require to continue at least until k > lambda
while (max(term) > eps | k < lambda) & k < maxnterms,
   k = k+1;
   i = find(term>eps | k < lambda);
   term = zeros(size(term));
   term(i) = lambda(i).^k*sqrt(k+c)/factorial(k).*exp(-lambda(i));
   Esqrt = Esqrt + term;
end
errest = term;
n_terms = k;
% This procedure fails (produces NaN's) for lambda > 95 apprx. So, we correct
% with an approximation (provisional, could be improved)
i = find(lambda>95);
Esqrt(i) = sqrt(lambda(i)-1/4+c);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
