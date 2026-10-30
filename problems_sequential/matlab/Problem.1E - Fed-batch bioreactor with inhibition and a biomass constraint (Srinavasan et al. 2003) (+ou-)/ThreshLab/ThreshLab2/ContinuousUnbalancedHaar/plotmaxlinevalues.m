
function plotmaxlinevalues(maxlines,l,w,options)

% plotmaxlinevalues - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     plotmaxlinevalues(maxlines,l,w)
%     plotmaxlinevalues(maxlines,l,w,options)
%   Inputs
%     maxlines   Matrix whose columns are plots of positions of local maxima as
%                function of scale
%     l          index to specific line in maxlines
%     w          matrix of continuous wavelet coefficients corresponding to
%                maxlines
%     options    Matlab plot options
%   Outputs
%     <none>
%   Description
%      Plots coefficient values along line (l) of local maxima in maxlines

if nargin < 4, options = 'b'; end
maxscales = find(maxlines(:,l));
maxx = maxlines(maxscales,l);
for s = 1:length(maxscales), values(s) = w(maxscales(s),maxx(s)); end
plot(maxscales,values,options)
