
function F = CDFmass0norm(v,p,mu1,sig1);

% usage: F = CDFmass0norm(v,p,mu1,sig1);
%
% evaluates cumulative distribution function of a random variable, say X, with
% point mass at zero and normal density away from zero.
%
% i.e., X|X=/=0 ~ n(mu1,sig1) and P(X=0) = 1-p
%
% input: v: value in which distribution function is evaluated
%        p: probability that the random variable is nonzero
%      mu1: E(X|X=/=0) (this expected value may be zero)
%     sig1: sqrt(V(X|X=/=0))

eps = 1.e-8;
q = 1-p;
if sig1 <eps,
   F = q.*(v >= 0) + p.*(v >= mu1);
else
   F = q.*(v >= 0) + p .* cumgauss(v-mu1,sig1);
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
