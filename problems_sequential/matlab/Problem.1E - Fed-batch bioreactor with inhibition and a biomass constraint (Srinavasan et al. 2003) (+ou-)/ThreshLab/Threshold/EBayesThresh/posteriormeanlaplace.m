
function postmean = posteriormeanlaplace(w,p,a,sigma);

% posteriormeanlaplace -- find posterior mean in
%                           zero-inflated Laplace + normal noise model
%   Usage
%     vhat = posteriormeanlaplace(w,p,a,sigma);
%   Inputs
%     w     observation
%     p     prior probability P(V =/= 0)
%     a     parameter of Laplace distribution (default 0.5)
%     sigma st.dev. of normal distribution (default 1)
%   Outputs
%     vhat  posterior mean, i.e., vhat satisfies: vhat = E(V|W=w)
%   Description
%      The model is W = V + U, with
%      p = P(V =/= 0)
%      V|V=/=0 ~ Laplace(a)
%      U ~ NID(0,sigma^2)


if nargin < 4, sigma = 1; end
if nargin < 3, a = 1; end

s = sign(w);
w = abs(w);
v0 = zeros(size(w));
v1 = w;

eps = 1.e-308;
% take eps such that 1/eps is still finite. It may happen that eps > 0, but
% 1/eps yields Inf.
ppost = posteriorplaplace(p,w,a,sigma);
oo = find(gauss(abs(w)+sigma^2*a,sigma) < eps);
ww = w;
ww(oo) = 0;
postmean = ww + sigma*a^2/2 * gauss(ww,sigma) .* ...
  (Rcumgauss(ww+sigma^2*a,sigma)./gauss(ww+sigma^2*a,sigma) - ...
    cumgauss(ww-sigma^2*a,sigma)./gauss(ww-sigma^2*a,sigma));
postmean = s.*postmean.*ppost;
postmean(oo) = s(oo).*w(oo);
