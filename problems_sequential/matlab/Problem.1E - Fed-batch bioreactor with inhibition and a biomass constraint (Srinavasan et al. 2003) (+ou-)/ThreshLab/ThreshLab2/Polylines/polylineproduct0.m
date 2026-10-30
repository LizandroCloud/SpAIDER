
function P = polylineproduct0(F,G,xF,xG);

% polylineproduct0 - ThreshLab2/Polylines --- inner product of two polylines in
% L_2; fast computation for single polylines
%   Usage
%     P = polylineproduct0(F,G,xF,xG);
%     P = polylineproduct0(F,G,x);
%     P = polylineproduct0(F,G);
%   Inputs
%     F     first polyline values in knots (possibly a matrix)
%     G     second polyline values in knots (must not be a matrix)
%     x     (common) x-values of knots of both polylines
%     xF    x-coordinates of knots of first polyline
%     xG    x-coordinates of knots of second polyline
%   Outputs
%     P     inner product of polylines
%   Description
%     computes innerproducts of the columns of F with the columns of G along
%     the grid points given in xF and xG
%     This version can deal with matrix-vector products and inner products, but
%     not with matrix-matrix products: G must not be a matrix
%   See also
%     help polylineproduct
%     help polylineproduct1
%     help polyline

rF = size(F,1);
rG = size(G,1);
if rF == 1, F = F'; rF = size(F,1); end
if rG == 1, G = G'; rG = size(G,1); end
if nargin < 3, xF = (1:rF); end
if nargin < 4,
   if rF == rG, xG = xF; else xG = (1:rG); end
end
if rF ~= length(xF), F = F'; rF = size(F,1); end
if rG ~= length(xG), G = G'; rG = size(G,1); end

if any(xF ~= xG),
   x = column(unique([row(xF) row(xG)]));
   F = polyline(x,xF,F);
   G = polyline(x,xG,G);
else
   x = xF;
end

cG = size(G,2); cF = size(F,2);
if cG>1, FG = F; F = G; G = FG; end


N = length(x);
Delta = column(x(2:N)-x(1:N-1));
Delta1 = [Delta; 0];
Delta2 = [0; Delta];

G1 = G.*Delta1;
G2 = G.*Delta2;

P = (F'*G1 + F'*G2)/3 + (F'*shift(G2) + F'*shift(G1,-1))/6;
% Zie schr. 15.45
