
function [even, odd] = splitevenodd(x,varargin)

%  splitevenodd -- partitions data into even and odd indexed ordered subsets
%   Usage
%     [even, odd] = splitevenodd(x);
%     [even, odd] = splitevenodd(N);
%   Inputs
%     x      data vector, grid location vector, ...
%     N      length of (undisclosed) data vector
%   Outputs
%     even   ordered indices of even subset
%     odd    ordered indices of odd subset
%   Description
%   Note
%     other input arguments are ignored
%   See also
%     help FWT_2GPpolintU
%     help FWT_2GUPpolint
%     (and other second generation wavelet transforms)

N = length(x);
if N == 1, N = x; end
even = (1:2:N);
odd = (2:2:N);
