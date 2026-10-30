
function plotmaxlines(maxlines,options)

% plotmaxlines - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     plotmaxlines(maxlines)
%     plotmaxlines(maxlines,options)
%   Inputs
%     maxlines   Matrix whose columns are plots of positions of local maxima as
%                function of scale
%     options    Matlab plot options
%   Outputs
%     <none>
%   Description
%      Plots all lines of local maxima in maxlines

if nargin < 2, options = 'b'; end
h = ishold;
for l = 1:size(maxlines,2),
   maxscales = find(maxlines(:,l));
   maxx = maxlines(maxscales,l);
   plot(maxx,maxscales,options)
   hold on
end
if ~h, hold off, end

