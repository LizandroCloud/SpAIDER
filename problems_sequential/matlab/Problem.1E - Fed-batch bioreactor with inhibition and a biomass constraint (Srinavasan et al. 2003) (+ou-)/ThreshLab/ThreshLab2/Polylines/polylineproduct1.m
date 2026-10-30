
function P = polylineproduct1(F,G,xF,xG);

% polylineproduct1 - ThreshLab2/Polylines --- inner product of two polylines in
% L_2; version for matrices F and G
%   Usage
%     P = polylineproduct1(F,G,xF,xG);
%     P = polylineproduct1(F,G,x);
%     P = polylineproduct1(F,G);
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
%     This version can deal with matrices F and G, but is supotimal wrt 
%     computation time for vectors
%   See also
%     help polylineproduct
%     help polylineproduct0
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

N = length(x);
Delta = sparse(diag(x(2:N)-x(1:N-1)));

F1 = F(1:N-1,:);
G1 = G(1:N-1,:);
F2 = F(2:N,:);
G2 = G(2:N,:);

P = (F1'*Delta*G1 + F2'*Delta*G2)/3 + (F1'*Delta*G2 + F2'*Delta*G1)/6;
% Zie schr. 15.45
