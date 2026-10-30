
function maxlines = removelowmaxlines(maxlines,z,threshold)

% removelowmaxlines - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     maxlines = removelowmaxlines(maxlines,z,threshold);
%   Inputs
%     maxlines   Matrix whose columns are plots of positions of local maxima as
%                function of scale
%     z          matrix of continuous wavelet transform coefficients,
%                corresponding to maxlines
%     threshold  (optional; default = 3) threshold of significance
%   Outputs
%     maxlines   Modified version of input matrix. Dimensions may have changed.
%   Description
%     Once lines of local maxima have been constructed, it should be checked
%     whether at least any of the maxima on a line could be significant (i.e.,
%     above a threshold of significance. Lines that do not reach this threshold
%     are removed.
%   See also
%      help connectmaxlines
%      help deconnectmaxlines
%      help mergemaxlines

if nargin<3,
   threshold = 3;
end

ll = 0;
maxllines = zeros(size(maxlines));

for l = 1:size(maxlines,2),
   maxscales = find(maxlines(:,l));
   % maxscales are levels (j) where a maximum of this line (l) is present
   maxx = maxlines(maxscales,l);
   clear maxz maxzz
   for j = 1:length(maxscales),
      maxz(j) = z(maxscales(j),maxx(j));
   end
   if max(abs(maxz)) > threshold,
      ll = ll+1;
      maxllines(:,ll) = maxlines(:,l);
      % i.e., only keep lines of local maxima if at least one maximum has an
      % absolute z-value equal or above threshold
   end
end

maxlines = maxllines(:,1:ll);
