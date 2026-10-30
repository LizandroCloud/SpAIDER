
function maxim = removeclosemaxima(maxim,z,M,m);

% removeclosemaxima - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     maxim = removeadjacentmaxima(maxim,z,M,m);
%   Inputs
%     maxim     binary pixel matrix indicating which coefficient in z is a
%               local maximum
%     z         matrix of continuous wavelet transform coefficients
%     M         number of dyadic scales in discretized continuous transform
%     m         number of additional intermediate scales (see help CHT)
%   Outputs
%     maxim     binary pixel matrix, with the same size as the input version
%   Description
%     It often happens that one feature + noise causes two local maxima in the
%     wavelet coefficients that are closer to each other than the resolution at
%     the level of the coefficients. The presence of two maxima closer to each
%     other than the resolution is due to noise, and therefore these maxima
%     should be reduced to a single one.
%   See also
%      help CHT
%      help findCWTmaxima

[dummy LR] = CHT(1,M,m);
LR = LR/2;
% LR contains for every scale the number of observations left and right
% involved in the computation of a balanced Haar wavelet coefficient at that
% scale.
Mm = length(LR); % Mm = size(z,1);
N = size(z,2);
for j = 1:Mm,
   scale = LR(j);
   for k = 1:N,
      if maxim(j,k),
         l = max(1,k-scale); r = min(N,k+scale);
         if abs(z(j,k)) < max(abs(z(j,l:r))), maxim(j,k) = 0; end
      end
   end
end
