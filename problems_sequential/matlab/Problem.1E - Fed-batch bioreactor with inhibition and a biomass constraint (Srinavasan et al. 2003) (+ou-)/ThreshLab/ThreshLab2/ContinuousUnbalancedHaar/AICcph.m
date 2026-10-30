
function aic = AICcph(y,cp,pdf);

% AICcph - ThreshLab2/ContinuousUnbalancedHaar Akaiki's Information Criterion
%                                            for multiple Change Point
%                                            selection using unbalanced,
%                                            continuous Haar wavelet analysis
%   Usage
%     aic = AICcph(y,cp,pdf);
%   Inputs
%     y          observations
%     cp         locations of change points
%     pdf        'normal' (default) or 'poisson'
%   Outputs
%     aic        AIC-value for given choice of cp

if nargin < 3, pdf = 'normal'; end

n = length(y);
cp = row(sort(cp));
if cp(1) > 1, cp = [1 cp]; end
ncp = length(cp);
if cp(ncp) < n+1, cp = [cp n+1]; ncp = ncp + 1; end
for b = 1:ncp-1,
   bp0 = cp(b);
   bp1 = cp(b+1);
   ff(bp0:bp1-1) = mean(y(bp0:bp1-1));
end

if strcmp(pdf,'poisson') | strcmp(pdf,'Poisson') | strcmp(pdf,'pois'),
   ff1 = ff.*double(ff>0) + double(ff==0);
   LL = sum(-ff + y.*log(ff1).*double(ff>0) - logfactorial(y));
   ff0 = find(ff==0);
   if max(y(ff0)) > 0, LL = Inf; end
else
   stdev = sqrt(sum((y-ff).^2)/n);
   LL = -n*log(stdev) - n*(log(2*pi)+1)/2;
end
aic = -2*LL + 2*(2*ncp-3); % 2*ncp - 3 is number of estimated parameters
                       % ncp - 1 is the trace of the regression matrix
   












% http://www-ssc.igpp.ucla.edu/personnel/russell/ESS265/Ch9/autoreg/node15.html#SECTION00043000000000000000
