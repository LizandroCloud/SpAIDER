
function [y, L] = dyadicextension(y,L);

% usage: [y L] = dyadicextension(y,L);
% or: y = dyadicextension(y);
%
% symmetric extension of a signal, followed by extension to have dyadic input
% vector.
% also make sure that the number of levels in the wavelet decomposition remains
% the same. Since the lowest level L has more scaling coefficients after
% extension of the signal, the value of L has to be adjusted (if this value
% is given)

N = length(y);
% first extend signal to avoid boundary effects
y = [y(N:-1:1) y y(N:-1:1)];
% next, extend signal to have dyadic number of data
N0 = N;
N = length(y);
J = ceil(log2(N));
yext = y(N-1:-1:2*N-2^J);
y = [y yext];

if nargin > 1,
   N = length(y)
   M = ceil(log2(N0))-L;
   L0 = L;
   J = log2(N);
   L = J-M;
end
