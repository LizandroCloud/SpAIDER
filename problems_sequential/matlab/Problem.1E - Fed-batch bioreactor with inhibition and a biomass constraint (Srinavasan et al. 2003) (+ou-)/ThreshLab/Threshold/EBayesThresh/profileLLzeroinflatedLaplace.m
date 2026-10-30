
function [LL a] = profileLLzeroinflatedLaplace(p,W);

% profileLLzeroinflatedLaplace -- profile log likelihood for zero inflated
% Laplace + standard normal observational noise
%   Usage
%      LL = profileLLzeroinflatedLaplace(p,W);
%  [LL a] = profileLLzeroinflatedLaplace(p,W);
%   Inputs
%      W   observations
%      p   parameter value
%   Outputs
%      LL  profile marginal log likelihood of p given W
%   Description
%      The model is W = V + Z, with
%      p = P(V =/= 0)
%      V|V=/=0 ~ Laplace(a)
%      Z ~ NID(0,1)

stdev = 1;

W2bar = mean(W.^2);
W2bar = max(W2bar,stdev^2);
ahat = 2*p/sqrt(W2bar-stdev^2);

LL = logmarginallikelihoodlaplace(p,ahat,1,W);
a = ahat;
