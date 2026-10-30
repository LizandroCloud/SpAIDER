
function s2 = varestGCV(x,soft)

% varestGCV: estimation of the variance, using GCV
%  Usage
%    s2 = varestGCV(x);
%    s2 = varestGCV(x,soft);
%  Inputs
%    x      input data
%    soft   (optional) binary value: 0 = hard thr.; 1 = soft thr.
%  Outputs
%    s2     variance estimation
%  Description
%    default: (soft = 1)
%    soft threshold shrinkage using GCV threshold, followed by variance
%    estimation based on residuals
%    alternative: (soft = 0)
%    hard thresholding using GCV threshold, followed by variance estimation
%    based on residuals
%  See also
%    help mingcvsoft
%    help GCVsoft

if nargin == 1, soft = 1; end
if soft,
   xt = mingcvsoft(x);
else
   xt = mingcvhard(x);
end
n0 = sum(abs(xt)<eps);
% s2 = sum((xt-x).^2)/length(x);
s2 = sum((xt-x).^2)/n0;
