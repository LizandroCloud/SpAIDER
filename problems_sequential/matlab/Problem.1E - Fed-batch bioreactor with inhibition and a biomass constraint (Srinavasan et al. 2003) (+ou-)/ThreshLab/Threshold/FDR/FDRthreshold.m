function thr = FDRthreshold(x,stdev,alfa)

% FDRthreshold -- finds false discovery rate threshold
%  Usage
%    thr = FDRthreshold(x,stdev,alfa);
%    thr = FDRthreshold(x,stdev);
%    thr = FDRthreshold(x);
%  Inputs
%    x      observations, assumed to be sparse and i.i.d. normal
%    stdev  standard deviation of the noise (or an estimate of it)
%    alfa   target false discovery rate (default 0.05)
%  Outputs
%    thr    threshold value
%  Description
%    finds thr such that FDR < alfa
%    if stdev is NaN, 0, Inf, or non-existent at function call, then it is
%    estimated using Median Absolute Deviation.

if nargin < 3, alfa = 0.05; end
if nargin < 2, stdev = NaN; end

if isnan(stdev) | stdev == 0 | isinf(stdev),
   stdev = MAD1(x)/invcumgauss(0.75);
end

if length(stdev) == 1, stdev = stdev*ones(size(x)); end

p = 2*Rcumgauss(abs(x)./stdev);

[psort i] = sort(p);
xsort = x(i);
n = length(p);
k = max(find(row(psort) - (1:n)/n*alfa < 0));
% plot(psort), hold on, plot((1:n)/n*alfa,'r'), hold off
if isempty(k), k = 1; end
thr = abs(xsort(k));
