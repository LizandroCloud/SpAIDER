
function Mn = initmom(x,n);

% initmom -- ThreshLab2/Irregular/1D/ -- 
%     initializes moments of scaling basis functions at finest scale
%  Usage
%    Mn = initmom(x,n);
%    Mn = initmom(x);
%  Inputs
%    x      grid-locations of observations (time, space)
%    n      vector of exponents (must not contain -1, default is 0)
%  Outputs
%    Mn     matrix of size length(x) times length(n)
%  Description
%    computes integral(phi(x)*x^n) with phi(x) a characteristic function on 
%    [a(i),a(i+1)], where a(1) = x(1) and a(i) = (x(i) + x(i+1))/2
%    This may serve as approximations for the finest scales scaling function
%    moments.
%  See also
%    help IRT_2GprefiltPkernelUpdate
%    help RWT_2GprefiltPpolintUpdate
%    help FWT_2GkernelBApyramid

if nargin < 2,
   n = 0;
end
x = row(x);
n = column(n);
p = length(n);
N = length(x);
Mn = zeros(p,N);
a(N+1) = x(N);
a(1) = x(1);
a(2:N) = (x(1:N-1)+x(2:N))/2;
a = repmat(a,size(n));
n = repmat(n,size(x));
Mn = (a(:,2:N+1).^(n+1) - a(:,1:N).^(n+1))./(n+1);
Mn = Mn';
