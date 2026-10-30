
function [muhat thr stdevest phat ahat] = ebayesthresh(x,varargin);

% ebayesthresh -- Empirical Bayesian Thresholding (Johnstone & Silverman, 2004)
%  Usage
%    [muhat thr phat ahat stdevest] = ebayesthresh(x,...)
%    [muhat thr phat ahat stdevest] = ebayesthresh(x,...
%        'prior',priormodel,...
%        'bayesfac',0_or_1,...
%        'stdev',stdev,...
%        'threshrule',threshrule,...
%        'isotone',0_or_1,...);
%        'a',a
%  Inputs
%    x         observations of X = mu + noise with mu sparse
%    'prior',priormodel
%              prior model for mu; 'Laplace','DE' or 'Cauchy'
%    'bayesfac',0_or_1
%              binary value indicating whether or not to use bayes factors
%    'stdev',stdev
%              (estimation of) standard deviation
%    'threshrule',threshrule
%              'postmedian': posterior median
%              'postmean'  : posterior mean
%    'isotone',0_or_1
%  Outputs
%    muhat     thresholded values of x, estimates for mu = E(X)
%    thr       the (implicit) threshold (if applicable)
%    stdevest  estimate of noise standard deviation
%    phat      estimate of p, degree of sparsity p = P(X =/= 0)
%    ahat      estimate of parameter in Laplace prior (if applicable)
%  NOTE
%    this routine is UNDER CONSTRUCTION, several options (such as isotone) have
%    not been implemented yet.

prior = 'laplace';
a = 1/2;
bayesfac = 0;
stdev = 1;
threshrule = 'median';
isotone = 0;

ia = 2;
while ia < nargin,
   switch varargin{ia-1},
   case {'prior'}
      prior = varargin{ia}; ia = ia+2;
   case {'bayesfac'},
      bayesfac = varargin{ia}; ia = ia+2;
   case {'stdev'},
      stdev = varargin{ia}; ia = ia+2;
   case {'threshrule'},
      threshrule = varargin{ia}; ia = ia+2;
   case {'isotone'},
      isotone = varargin{ia}; ia = ia+2;
   case {'a'},
      a = varargin{ia}; ia = ia+2;
   otherwise
      ia = ia+1;
   end
end

if isempty(stdev) | isnan(stdev) | isinf(stdev),
   stdev = MAD(x)/invcumgauss(3/4);
end
stdevest = stdev;

if length(stdev) == 1,
   stdev = stdev*ones(size(x));
end
if size(stdev) == size(x),
   x = x./stdev;
else
   error('dimension mismatch observations (x) and st.dev.')
end
  

switch prior(1:2),
case {'La','la','DE','de'}
   [phat ahat] = estimateparzeroinflatedLaplace(x,a);
   switch threshrule
   case {'posmedian','postmedian','posteriormedian','median','med', 'Median'},
      muhat = posteriormedianlaplace(x,phat,ahat,stdev);
   case {'posmean','postmean','posteriormean','mean','Mean'},
      muhat = posteriormeanlaplace(x,phat,ahat,stdev);
% case {'hard', 'Hard' , 'HT'}
% case {'soft', 'Soft', 'ST'}
% 
   case 'none',
      muhat = NaN;
   otherwise
      muhat = NaN
   end
case {'Ca','ca','JS'}
   phat = estimateparzeroinflatedCauchylike(x);
   phat = 1;
end
   

muhat = muhat.*stdev;
