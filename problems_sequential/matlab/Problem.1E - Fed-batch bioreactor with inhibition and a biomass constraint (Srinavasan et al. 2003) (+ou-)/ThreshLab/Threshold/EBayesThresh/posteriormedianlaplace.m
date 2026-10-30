
function postmedian = posteriormedianlaplace(w,p,a,sigma);

% posteriormedianlaplace -- find posterior median in
%                           zero-inflated Laplace + normal noise model
%   Usage
%     vhat = posteriormedianlaplace(w,p,a,sigma);
%   Inputs
%     w     observation
%     p     prior probability P(V =/= 0)
%     a     parameter of Laplace distribution (default 0.5)
%     sigma st.dev. of normal distribution (default 1)
%   Outputs
%     vhat  posterior median, i.e., vhat satisfies
%           F_V(vhat|W=w) = P(V < vhat|W=w) = 1/2
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

eps = 1.e-6;
ppost = posteriorplaplace(p,w,a,sigma);
p0 = ppost.*posteriorcdflaplace(v0,w,a,sigma) + 1-ppost;
p1 = ppost.*posteriorcdflaplace(v1,w,a,sigma) + 1-ppost;
v1(p0> 0.5) = 0;
notconverged = find(abs(v1-v0)>eps);
while ~isempty(notconverged),
   vv = (v0(notconverged)+v1(notconverged))/2;
   ww = w(notconverged);
   pp = ppost(notconverged).*posteriorcdflaplace(vv,ww,a,sigma) + ...
           1-ppost(notconverged);
   v1(notconverged(pp>=0.5)) = vv(pp>=0.5);
   v0(notconverged(pp<0.5)) = vv(pp<0.5);
   notconverged = find(abs(v1-v0)>eps);
end
postmedian = s.*(v0+v1)/2; 
