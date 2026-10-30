
function F = posteriorcdflaplace(v,w,a,sigma);

% posteriorcdflaplace -- find posterior cumulative distribution function in 
%                        zero-inflated Laplace + normal noise model
%   Usage
%     F = posteriorcdflaplace(v,w,a,sigma);
%   Inputs
%     p     prior probability P(V =/= 0)
%     w     observation
%     a     parameter of Laplace distribution (default 0.5)
%     sigma st.dev. of normal distribution (default 1)
%   Outputs
%     F     posterior cumulative distribution function in v, i.e.,
%           F_V(v|W=w) = P(V < v|W=w)
%   Description
%      The model is W = V + U, with
%      p = P(V =/= 0)
%      V|V=/=0 ~ Laplace(a)
%      U ~ NID(0,sigma^2)
%   Note
%      The expression used in this routine holds for v*w > 0 only

if nargin < 4, sigma = 1; end
if nargin < 3, a = 1; end

if size(w) ~= size(v),
   nw = length(w);
   nv = length(v);
   v = repmat(column(v),1,nw);
   w = repmat(row(w),nv,1);
end

eps = 1.e-308;
% take eps such that 1/eps is still finite. It may happen that eps > 0, but
% 1/eps yields Inf.
oo = find(gauss(w-sigma^2*a,sigma)<eps);
ww = w;
ww(oo) = 0;

F = 1 - (a/2)*gauss(ww,sigma)./laplaceplusnormalpdf(ww,a,sigma).*...
         cumgauss(ww-v-sigma^2*a,sigma)./gauss(ww-sigma^2*a,sigma);
F2 = 1 - cumgauss(w-v-sigma^2*a,sigma);
F(oo) = F2(oo);
