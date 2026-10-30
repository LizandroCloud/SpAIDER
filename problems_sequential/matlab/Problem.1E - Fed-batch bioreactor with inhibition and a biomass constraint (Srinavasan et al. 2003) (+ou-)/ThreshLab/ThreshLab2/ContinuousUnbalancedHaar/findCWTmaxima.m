
function [maxim ws] = findCWTmaxima(w,s);

% findCWTmaxima - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     [maxim ws] = findCWTmaxima(w,s);
%     maxim = findCWTmaxima(w,s);
%     maxim = findCWTmaxima(w,S);
%   Inputs
%     w       continuous wavelet transform coefficients
%     s       smoothing filter (default is 1, i.e., no smoothing)
%     S       smoothing matrix (dimensions should fit with w; default = eye)
%   Outputs
%     maxim   binary image with local maxima at each scale
%     ws      smoothed coefficients
%   Description
%     finds local maxima within each scale
%     if w is noisy, then preprocessing by filtering with s or S is possible
%   See also
%     help localmaxima
%     help connectmaxlines

if nargin < 2, s = 1; end
maxim = zeros(size(w));
M = size(w,1);
N = size(w,2);
ws = w;
for j = 1:M,
   wj = w(j,1:N);
   wjs = wj;
   if min(size(s)) == 1 & length(s) > 1,
      [m mp] = max(s);
      wjs = conv(s,wj);
      wjs = wjs(mp:mp+N-1);
   elseif size(s,1) == size(s,2) & size(s,1) == length(wj),
      wjs = s*column(wj);
      if size(wjs) ~= size(wj), wjs = wjs'; end
   end
   k = localmaxima(abs(wjs));
   % first find maxima in smoothed sequence of coefficients, then find
   % precise location of maxima in unsmoothed sequence
   kk = [1 row(k) N];
   kk = floor((kk(1:length(kk)-1) + kk(2:length(kk)))/2);
   for l = 1:length(k),
      % [wk k(l)] = max(abs(w(j,kk(l):kk(l+1))));
      % working with abs. values induces false local values at boundaries.
      [wk k(l)] = max(w(j,kk(l):kk(l+1))*sign(w(j,k(l))));
      k(l) = k(l) - 1 + kk(l);
   end
   maxim(j,k) = 256;
   ws(j,1:N) = wjs;
end
