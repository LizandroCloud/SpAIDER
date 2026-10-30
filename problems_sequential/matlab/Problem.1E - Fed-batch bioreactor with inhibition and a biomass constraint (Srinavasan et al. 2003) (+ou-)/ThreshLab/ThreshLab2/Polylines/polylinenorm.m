
function P = polylinenorm(f,x);

% polylinenorm - ThreshLab2/Polylines --- L_2-norm of a polyline defined by knots
%   Usage
%     P = polylinenorm(f,x);
%     P = polylinenorm(f);
%   Inputs
%     F     polyline values in knots
%     x     x-values of knots
%   Outputs
%     P     L_2-norm of a polyline
%   Description
%     computes sqrt(polylineproduct(f,f,x));
%   See also
%     help polylineproduct
%     help polyline

[r c] = size(f);
if r == 1, f = f'; end
[r c] = size(f);

if nargin < 2, x = (1:r); end

P = sqrt(polylineproduct(f,f,x));
