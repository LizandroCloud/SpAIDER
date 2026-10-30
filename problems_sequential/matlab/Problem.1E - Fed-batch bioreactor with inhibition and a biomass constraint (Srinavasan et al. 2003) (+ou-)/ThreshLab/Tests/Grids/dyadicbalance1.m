
function [beta, j, k] = dyadicbalance1(x);

% homogeneity1 - ThreshLab/Tests/Grids/ --- Dyadic balance for 1D grids
%   Usage
%     [beta, j, k] = dyadicbalance1(x);
%     beta         = dyadicbalance1(x);
%   Inputs
%     x            data points
%   Outpouts
%     beta         Dyadic balance constant (see definition below)
%     j, k         scale and location of most unbalanced subdivision
%   Description
%     This function computes the dyadic balance of a grid x, defined as follows.
%     Let Delta{J} = x(2:N) - x(1:N-1), 
%     and Delta{j} = Delta{j+1}(even) + Delta{j+1}(odd)
%     then
%     beta(j) = min(min(Delta{j+1}(even),Delta{j+1}(odd))./Delta{j})
%     beta = min(beta)
%     It holds: beta >= 1/(1+gamma), with gamma homogeneity constant (see
%     help homogeneity1)
%     maximum value for beta is 1/2, this is a equidistant grid
%     minimum value for beta is 0, such a grid has homogeneity constant Inf
%   See also
%     help homogeneity1
%   Reference
%     Ingrid Daubechies, Igor Guskov, and Wim Sweldens, Regularity of Irregular
%     Subdivision, Constructive Approximation, 15, pages 381-426, 1999.


x = sort(x);
N = length(x);
if length(x) < 3,
   beta = NaN;
else
   beta = 1/2;
   deltafine = x(2:N) - x(1:N-1);
   Nf = length(deltafine);
   J = floor(log2(Nf));
   j = 0; k = 0;
   for jj = J-1:-1:1,
      Nc = floor(Nf/2); even = (1:2:2*Nc); odd = (2:2:2*Nc);
      deltacoarse = deltafine(even) + deltafine(odd);
      deltamin = min(deltafine(even),deltafine(odd));
      [betaj kj] = min(deltamin./deltacoarse);
      deltafine = deltacoarse; Nf = Nc;
      if betaj < beta,
         beta = betaj;
         j = jj;
         k = kj(1);
      end
   end
end
