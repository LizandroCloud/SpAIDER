
function P = polylineproduct(F,G,xF,xG);

% polylineproduct - ThreshLab2/Polylines --- inner product of two polylines in
% L_2
%   Usage
%     P = polylineproduct(F,G,xF,xG);
%     P = polylineproduct(F,G,x);
%     P = polylineproduct(F,G);
%   Inputs
%     F     first polyline values in knots
%     G     second polyline values in knots
%     x     (common) x-values of knots of both polylines
%     xF    x-coordinates of knots of first polyline
%     xG    x-coordinates of knots of second polyline
%   Outputs
%     P     inner product of polylines
%   Description
%     computes innerproducts of the columns of F with the columns of G along
%     the grid points given in xF and xG
%     This routine chooses between polylineproduct0 and polylineproduct1
%   See also
%     help polylineproduct0
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

mergeFG = 0;
if length(xF) ~= length(xG), mergeFG = 1; elseif any(xF ~= xG), mergeFG = 1; end
if mergeFG,
   x = column(unique([row(xF) row(xG)]));
   F = polyline(x,xF,F);
   G = polyline(x,xG,G);
else
   x = xF;
end
cG = size(G,2); cF = size(F,2);
if cG>1 & cF>1, 
   P = polylineproduct1(F,G,x,x);
else
   if cG>1, FG = F; F = G; G = FG; end
   P = polylineproduct0(F,G,x,x);
end
