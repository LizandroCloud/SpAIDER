
function LL = logmarginallikelihoodlaplace(pp,aa,sigma,w);

% logmarginallikelihoodlaplace -- marginal log-likelihood for zero inflated
% Laplace + standard normal observational noise
%   Usage
%      LL = logmarginallikelihoodlaplace(pp,aa,sigma,W);
%      LL = logmarginallikelihoodlaplace(pa,sigma,W);
%   Inputs
%      W     observations
%      pp    sparsity parameter value
%      aa    parameter Laplace density
%      pa    [pp aa] or [pp; aa]
%      sigma value for st.dev. of noise
%   Outputs
%      LL    marginal log likelihood of p given W
%   Description
%      The model is W = V + U, with
%      p = P(V =/= 0)
%      V|V=/=0 ~ Laplace(a)
%      U ~ NID(0,sigma^2)

if nargin < 4,
   w = sigma;
   sigma = aa;
   pa = pp;
   [s1 s2] = size(pa);
   if s2 == 2,
      pp = pa(:,1); aa = pa(:,2);
   elseif s1 == 2,
      pp = pa(1,:); aa = pa(:,2);
   else
      error('inconsistent input variables')
   end
end
na = length(aa); np = length(pp);
if na < np,
   aa = aa(1)*ones(size(pp));
elseif np < na,
   pp = pp(1)*ones(size(aa));
end
np = length(pp); 
LL = zeros(size(pp));
for i = 1:np, 
   p = pp(i);
   a = aa(i);
   LL(i) = sum(log(p*laplaceplusnormalpdf(w,a,sigma) + (1-p)*gauss(w,sigma)));
end
