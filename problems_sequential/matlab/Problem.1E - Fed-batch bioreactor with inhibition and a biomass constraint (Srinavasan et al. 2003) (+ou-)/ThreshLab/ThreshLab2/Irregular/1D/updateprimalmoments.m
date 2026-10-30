
function U = updateprimalmoments(xodd,yodd,xeven,Modd,Meven);

% updateprimalmoments - ThreshLab2/Irregular/1D/ --- update or primal lifting
% step, designed for vanishing moments of primal (=synthesis) wavelet basis
%   Usage
%      yeven = updateprimalmoments(xodd,yodd,xeven,Modd,Meven);
%      U     = updateprimalmoments(xodd,I,xeven,Modd,Meven);
%   Inputs
%      xodd    x-coordinates of odd observations/odd fine scale coefficients
%      yodd    y-values in odd indexed locations
%      I       matrix (e.g., identity matrix of size = length(xodd) )
%      xeven   ``even'' indexed locations, i.e., locations that proceed to
%              coarse scale
%      Modd    matrix whose p-th column is (p-1)st moment of the current odd
%              scaling functions (first column thus contains integrals of
%              scaling functions). These moments are to be annihilated by
%              this routine, thereby turning the odd scaling functions into 
%              wavelet functions.
%      Meven   matrix whose p-th column is (p-1)st moment of the even scaling
%              functions (These scaling functions and their moments are left
%              untouched by this update operation)
%   Outputs
%      yeven  values to be added to in even indexed coefficients
%      U      for I = eye(length(xeven)), U is the update operation matrix
%   Description
%      This routine computes the update that annihilates the moments of the
%      primal (=synthesis) wavelet (= detail) basis functions in a lifting 
%      scheme. The number of moments that are annihilated (i.e., the number
%      of primal vanishing moments) equals the number of columns in Mfine and
%      Mcoarse. The columns of Mfine and Mcoarse provide sets of linear 
%      equations for the columns of U, hence, the number of nonzero elements in
%      each column of U is bounded by the number of columns in Mfine and
%      Mcoarse
%      The moment anihilation equation is
%
%         0 = Modd(:,p) - U'*Meven(:,p)
%
%      for p = 1,...,n  where n is the number of moments to be annihilated (the
%      number of vanishing moments)
%      This has to be interpreted as a set of n equations for every row (o) of
%      Modd. The unknowns of each set are the nonzero elements in the o-th
%      column of U
%
%   See Also
%      help polynomialprediction

n = size(Modd,2);
if size(Meven,2) ~= n, error('dimension mismatch in Modd and Meven'), end

Utranspose = 0;
if size(xeven,1) > 1, xeven = xeven'; Utranspose = 1; end
if size(xodd,1) > 1, xodd = xodd'; end

No = length(xodd); Ne = length(xeven);
if size(Meven,1) ~= Ne, 
   error('dimension mismatch in Meven and xeven')
end
if size(Modd,1) ~= No,
   error('dimension mismatch in Modd and xodd')
end



L = ceil(n/2); R = n-L;
% if n is odd, then one more even on the left is involved in the update then on
% the right (cfr. polynomialprediction)

U = zeros(Ne,No);
for o = 1:No,
   xo = xodd(o);
   e = max(find(xeven<=xo));
   e = (e-L+1:e+R); e = e((e>0)&(e<Ne+1));
   xe = xeven(e); kk = length(e);
   U(e,o) = (Meven(e,1:kk)')\(Modd(o,1:kk)');
end
U = sparse(U);
if size(yodd,1) ~= size(U,2), yodd = yodd'; end
U = U*yodd;
