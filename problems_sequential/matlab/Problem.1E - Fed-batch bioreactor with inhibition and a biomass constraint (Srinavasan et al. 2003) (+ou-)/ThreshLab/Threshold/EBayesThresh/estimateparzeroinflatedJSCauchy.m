
function phat = estimateparzeroinflatedJSCauchy(w);

% estimateparzeroinflatedJSCauchy -- estimate parameters p in model
% W = V + Z, with V zero inflated JSCauchy and Z standard normal noise
%  Usage
%     phat = estimateparzeroinflatedJSCauchy(w);
%  Inputs
%     w      observations
%  Outputs
%     phat   estimate of p (degree of sparsity)
% Description
%     The model is W = V + Z, with
%     p = P(V =/= 0)
%     V|V=/=0 ~ JSCauchy(1)
%     Z ~ NID(0,1)
%   Reference
%     The Cauchy-like distribution JScauchy comes from:
%     I. Johnstone and B. Silverman, Needles and Straw in Haystacks, Annals of
%                                    Statistics, 2004.

eps = 0.001; %i.e. 10 function evaluations
p0 = 0;
p1 = 1;
interval = [p0 p1];
L0prime = marginalscoreJScauchy(p0,1,w);
L1prime = marginalscoreJScauchy(p1,1,w);
if (L0prime*L1prime) > 0
   L0 = logmarginallikelihoodJScauchy(p0,1,w);
   L1 = logmarginallikelihoodJScauchy(p1,1,w);
   p = double(L0 < L1);
else
  while (p1-p0) > eps
    p = (p0+p1)/2;
    Lprime = marginalscoreJScauchy(p,1,w);
    if ((Lprime*L0prime) > 0)
        p0 = p;
        L0prime = Lprime;
    else
        p1 = p;
        L1prime = Lprime;
    end
  end
p = (p0+p1)/2;
end

phat = p;
