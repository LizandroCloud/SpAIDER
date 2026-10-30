
function y = interpolationneville(x,x0,y0);

% interpolationneville - ThreshLab2/Irregular/1D/ --- Neville's polynomial
% interpolation algorithm
%   Usage
%      y = interpolationneville(x,f,t);
%   Inputs
%      x  points in which interpolating polynomial is evaluated
%      x0 x-coordinates of interpolation points
%      y0 y-coordinates of interpolation points
%   Outputs
%      y  
%   Description
%      y0 may be a matrix (check sizes for consistency with x0)
%      x may be a vector
%      If x is a scalar, then y is a column vector if x0 is a row vector and
%      vice versa
%      If x is a row vector, then the number of columns in y equals the number
%      of columns in x
%      If x is a column vector, then the number of rows in y equals the number
%      of rows in x
%      The other dimension in y is determined by the dimension of y0 not
%      corresponding to x0
%   See also
%      help polynomialprediction

[m n] = size(x0);
x0transpose = 0;
if m > n, x0 = x0'; y0 = y0'; [m n] = size(x0); x0transpose = 1; end
if m > 1, error('input argument x0 must be a scalar or a 1D-vector'), end

[r s] = size(y0);
if ~(n == s),
   if ~(n == r),
      error('dimension mismatch: inputs x0 and y0')
   else
      y0 = y0'; [r s] = size(y0);
   end
end
% x0 is reshaped as a row vector; y0 is reshaped accordingly, such that all
% *rows* of y0 match the dimension of x0.
% x is also reshaped as a row vector

[M N] = size(x);
finaltranspose = 0;
if M > N, finaltranspose = 1; x = x'; [M N] = size(x); end
if M > 1, error('input argument x must be a scalar or a 1D-vector'), end
if N == 1, finaltranspose = x0transpose; end

for j = 1:n, y00{j} = repmat(y0(:,j),1,N); end
xx = repmat(x,r,1);
for i = 1:n-1
  for j = 1:n-i
     y00{j} = ...
        ( (xx-x0(j+i)).*y00{j} - (xx-x0(j)).*y00{j+1} ) / ( x0(j) - x0(j+i) );
  end
end
if isempty(y00)
  y = 0;
else
  y = y00{1};
end

if finaltranspose, y = y'; end
