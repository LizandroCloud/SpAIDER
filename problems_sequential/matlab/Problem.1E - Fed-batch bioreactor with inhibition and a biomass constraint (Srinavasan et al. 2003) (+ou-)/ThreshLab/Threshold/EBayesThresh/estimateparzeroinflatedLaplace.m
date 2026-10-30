
function [phat ahat] = estimateparzeroinflatedLaplace(w,a);

% estimateparzeroinflatedLaplace -- estimate parameters a and p in model
% W = V + Z, with V zero inflated Laplace and Z standard normal noise
%  Usage
%     [phat ahat] = estimateparzeroinflatedLaplace(w);
%     phat = estimateparzeroinflatedLaplace(w,a);
%  Inputs
%     w      observations
%  Outputs
%     ahat   estimate of a
%     phat   estimate of p
% Description
%     The model is W = V + Z, with
%     p = P(V =/= 0)
%     V|V=/=0 ~ Laplace(a)
%     Z ~ NID(0,1)


if nargin == 1, a = -1; end

eps = 1.e-4;
n = length(w);
stdev = 1;

w2bar = mean(w.^2);
w2bar = max(w2bar,stdev^2);
univthr = sqrt(2*log(n))*stdev;

[wthrGCV thrGCV] = mingcvhard(w,'maxthr',univthr);
[wthrSURE thrSURE] = minsurehard(w);

pmin = sum(abs(w)>univthr)/n;
% pmin is truely a lower bound for p, beacuse univthr is a very conservative
% threshold AND sum(abs(w)>univthr) only counts v whose observation is above
% the threshold, whereas many nonzero v have an observation below the
% threshold



if isempty(a) | isnan(a) | isinf(a) | a < 0,

   % first find initial guess, using profile likelihood minimization

   amax = sqrt(2); % then EV^2 = EZ^2 = 1
   pmax = pmin/(2*(1-laplaceplusnormalcdf(univthr,amax,1)));
   pmax = min(pmax,1);
   % profileLLzeroinflatedLaplace is not concave, so we first find the minimum
   pmax2 = minfunc1Dfibon('profileLLzeroinflatedLaplace','lo',pmin,'up',pmax,w);
   if pmax2 > pmin, pmax = pmax2; end

   phat0 = maxfunc1Dfibon('profileLLzeroinflatedLaplace','lo',pmin,'up',pmax,w);

   ahat0 = 2*phat0/sqrt(w2bar-stdev^2);
   amin = 2*pmin/sqrt(w2bar-stdev^2);
   amax = 2*pmax/sqrt(w2bar-stdev^2);

   amin = maxfunc1Dfibon('logmarginallikelihoodlaplace','argpos',2,...
         'lo',amin,'hi',amax,pmin,1,w);
   amax = maxfunc1Dfibon('logmarginallikelihoodlaplace','argpos',2,...
         'lo',amin,'hi',amax,pmax,1,w);
   pahat = optfunc2Dinvecfibon('logmarginallikelihoodlaplace',...
        'starting',[phat0 ahat0],'lo',[pmin amin],'hi',[pmax amax],'max',1,w);
   phat = pahat(1);
   ahat = pahat(2);

else

   ahat = a;

   % phat0 is value of p such that observed mean(W^2) equal stdev^2+2/a^2
   phat0 = (w2bar-stdev^2)*ahat^2/2;
   phat0 = min(phat0,1);
   % phat1 is value of p such that with given a the observed number of values
   % above universal threshold can be expected.
   phat1 = pmin/(2*(1-laplaceplusnormalcdf(univthr,a,1)));
   % pmin2 is the value of p such that postmedian threshold = univ. thr.
   pmin2 = a*cumgauss(univthr-a)/gauss(univthr-a) -  ...
           laplaceplusnormalpdf(univthr, a)/gauss(univthr) + 1;
   pmin2 = 1/pmin2;
   pmin = max(pmin,pmin2);
   pmax = min(1,2*max([phat0 phat1 pmin]));

   phat = maxfunc1Dfibon('logmarginallikelihoodlaplace',...
                         'lo',pmin,'hi',pmax,ahat,1,w);
end
