function Nj = numbercoefateachscale(N,M);

% numbercoefateachscale - ThreshLab2/Irregular/1D/ --- finds number of
% coefficients at each scale in a critically downsampled lifted wavelet
% transform
%   Usage
%      Nj = numbercoefateachscale(N,M);
%      Nj = numbercoefateachscale(x,M);
%   Inputs
%      N      number of grid points/observations
%      x      grid, observational vector
%      M      number of levels in decomposition
%   Outputs
%      Nj     vector of number of elements at each scale
%             Nj(1) is the number of scaling coefficients
%   See also
%      help splitevenodd

if length(N)>1, N = length(N); end

even = (1:N); Ne = N;
Nj = zeros(1,M+1);
for j = 1:M,
   even = splitevenodd(even); Ne = length(even); Nj(M+2-j) = N-Ne; N = Ne;
end
Nj(1) = Ne;
