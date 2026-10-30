
function dLLdp = marginalscoreJScauchy(p,sigma,w);

% marginalscoreJScauchy -- derivative of marginal log-likelihood for zero 
% inflated JScauchy + standard normal observational noise
%   Usage
%      dLLdp = marginalscoreJScauchy(p,sigma,W);
%      dLLdp = marginalscoreJScauchy(p,W);
%   Inputs
%      W     observations
%      p     sparsity parameter value
%      sigma value for st.dev. of noise (default 1)
%   Outputs
%      dLLdp marginal score, i.e., derivative of marginal log likelihood of p 
%            given W
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

% beta = ( JScauchyplusnormalpdf(w) - gauss(w) )./gauss(w);
ww = w+(w==0);
% the following is to avoid overflow, but introduces possible inaccuracies for
% very small p
ww = min(ww,10*sigma);
ww = max(ww,-10*sigma);
beta0 = sigma/2-1;
beta1 = sigma^2*(exp(ww.^2/(2*sigma^2)) - 1)./ww.^2 - 1;
beta = beta0.*(w == 0) + beta1.*(w ~= 0);
dLLdp = sum(beta./(1 + p*beta));
% dLLdp = sum(1./(p + 1./beta));

