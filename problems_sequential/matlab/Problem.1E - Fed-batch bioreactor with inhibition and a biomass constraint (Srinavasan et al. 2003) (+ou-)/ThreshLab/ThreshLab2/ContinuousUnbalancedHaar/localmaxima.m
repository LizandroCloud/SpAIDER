
function x = localmaxima(f);

% localmaxima - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     x = localmaxima(f);
%   Inputs
%     f     vector of function values
%   Outputs
%     x     binary vector, 1 corresponds to local maximum in f
%   Description
%     finds all local maxima in f
%   See also
%      help connectmaxlines


n = length(f);
ff = [min(f)-1 row(f) min(f)-1];
x = find( (ff(2:n+1) >= ff(1:n)) & (ff(2:n+1) >= ff(3:n+2)));
