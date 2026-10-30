
function ppost = posteriorplaplace(p,w,a,sigma);

% posteriorplaplace -- find posterior probability P(V =/= 0|W) in zero-inflated
%                      Laplace + normal noise model
%   Usage
%     ppost = posteriorplaplace(p,w,a,sigma);
%   Inputs
%     p     prior probability P(V =/= 0)
%     w     observation
%     a     parameter of Laplace distribution (default 0.5)
%     sigma st.dev. of normal distribution (default 1)
%   Outputs
%     ppost posterior probability P(V =/= 0|W)
%   Description
%      The model is W = V + U, with
%      p = P(V =/= 0)
%      V|V=/=0 ~ Laplace(a)
%      U ~ NID(0,sigma^2)


if nargin < 4, sigma = 1; end
if nargin < 3, a = 1/2; end

if size(w) ~= size(p),
   nw = length(w);
   np = length(p);
   p = repmat(column(p),1,nw);
   w = repmat(row(w),np,1);
end

ppost = p.*laplaceplusnormalpdf(w,a,sigma)./...
       (p.*laplaceplusnormalpdf(w,a,sigma) + (1-p).*gauss(w,sigma));
