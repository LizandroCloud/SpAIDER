function doublerows = finddoublerowsinmatrix(A,eps)

% doublerows - ThreshLab2/ --- find double rows in matrix
%   Usage
%    doublerows = finddoublerowsinmtarix(A,eps);
%   Inputs
%    A    matrix
%    eps  working precision: if |a-b| < eps, then a equals b in this routine
%         (default 1.e-10)
%   Outputs
%    doublerows  cell array of double rows

if nargin < 2, eps = 1.e-10; end

[nrows ncols] = size(A);
if nrows == 1, A = A'; nrows = ncols; ncols = 1; end

firstcolumn = A(1:nrows,1);
[firstcolumnsorted order] = sort(firstcolumn);
A = A(order,1:ncols);
sentinel = Inf*ones(1,ncols);
A = [A;sentinel];
doublerow = zeros(nrows,1);

row = 1;
while row < nrows,
   nextrow = row+1;
   while max(abs(A(row,1:ncols)-A(nextrow,1:ncols))) < eps,
       doublerow(nextrow) = 1;
       nextrow = nextrow + 1;
   end
   row = nextrow;
end

rowswithdoubles = find(doublerow(2:nrows)-doublerow(1:nrows-1)==1);
doublerows = {};
for r = 1:length(rowswithdoubles),
   row = rowswithdoubles(r);
   lastrow = row+1;
   while doublerow(lastrow+1), lastrow = lastrow + 1; end
   doublerows{r} = order(row:lastrow);
end
