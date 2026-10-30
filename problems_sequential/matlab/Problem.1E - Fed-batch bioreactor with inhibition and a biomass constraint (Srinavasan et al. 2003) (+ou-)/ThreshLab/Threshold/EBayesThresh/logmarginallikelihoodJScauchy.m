
function LL = logmarginallikelihoodJScauchy(pp,sigma,w);

% logmarginallikelihoodJScauchy -- marginal log-likelihood for zero inflated
% JScauchy + standard normal observational noise
%   Usage
%      LL = logmarginallikelihoodJScauchy(pp,sigma,W);
%      LL = logmarginallikelihoodJScauchy(pp,W);
%   Inputs
%      W     observations
%      pp    sparsity parameter value
%      sigma value for st.dev. of noise (default 1)
%   Outputs
%      LL    marginal log likelihood of p given W
%   Description
%      The model is W = V + U, with
%      p = P(V =/= 0)
%      V|V=/=0 ~ JScauchy(sigma)
%      U ~ NID(0,sigma^2)
%   Reference
%     The Cauchy-like distribution JScauchy comes from:
%     I. Johnstone and B. Silverman, Needles and Straw in Haystacks, Annals of
%                                    Statistics, 2004.

if nargin < 3,
   w = sigma;
   sigma = 1;
end

qq = 1-pp;

np = length(pp);
LL = zeros(size(pp));
for i = 1:np,
   p = pp(i); q= 1-p;
   LL(i) = sum(log(p*JScauchyplusnormalpdf(w,sigma) + q*gauss(w,sigma)));
end
