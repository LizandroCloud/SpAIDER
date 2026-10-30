
function [Ht Ht0] = findbiorthfilter(P,M);

% Ht = findbiorthfilter(P,M);
% finds a filter Ht, as sparse as possible, such that
% 1) Ht*P = I
% 2) M'*(I-P*Ht) = 0
%
% This can be seen as solving a system with boundary conditions. (1) is the
% non-homogeneous system equation. We find a sparse particular solution for (1)
% and a sparse complementary solution of the homogeneous equation Ht*P = 0
%
% We first find sparse solutions for Ht*P = I and Ht0*P = 0. We then find a
% sparse matrix A such that Ht+A*Ht0 satisfies M'*(I-P*Ht) = 0
%
% The procedure uses the fact that Ht*P = I has a solution Ht with bandwidth
% smaller than Ne = size(P,2). This in its turn assumes that the bandwidth of P
% is smaller than Ne. (Otherwise the existence of a solution Ht with bandwidth
% below Ne is not guaranteed)
%
% The homogeneous solution has the same size as the particular solution. In
% general, the rank of the homogeneous solution is often different from that of
% the particular solution. A necessary condition to be sure that Ht0 generates
% all possible homogeneous solutions: its size must be above the dimension of 
% the general homogeneous solution, so dim null space = n-Ne must be smaller
% than Ne, i.e., Ne >= n/2
% The reason for working with equally sized homogeneous and particular
% solutions is that it facilitates imposing the sparsity, using similar
% procedures for finding both, and using squaresparsesolveAXBisC

[N Ne] = size(P);
p = 0;
if nargin > 1,
   [p Nm] = size(M);
end

% First we find the midpoint of each row in Ht, i.e., the column around which
% the nonzeros will be concentrated. This corresponds to the central row in the
% corresponding column of P
nonzerosP = (abs(P)>eps);
bandwidth = zeros(1,Ne);
periodic = 0;
for j = 1:Ne,
   nonzeroj = find(nonzerosP(1:N,j));
   nz1 = min(nonzeroj); nz2 = max(nonzeroj);
   bandwidth(j) = nz2-nz1;
   midpoint(j) = floor((nz2+nz1)/2);
   if bandwidth(j) == N-1,
      periodic = 1;
      shiftednonzerosj = shift(nonzerosP(1:N,j),floor(N/2));
      nonzeroj = find(shiftednonzerosj);
      nz1 = min(nonzeroj); nz2 = max(nonzeroj);
      midpoint(j) = mod(floor((nz2+nz1)/2) - floor(N/2) -1,N)+1;
      % midpoint(j) = mod(midpoint(j)+floor(N/2)-1,N)+1;
   end
end

% Now, find sparse homogeneous and particular solutions Ht0 and Ht
Ht = zeros(size(P))';
Ht0 = zeros(size(Ht));
% rr will refer to rows of P that will be used for filling in nonzero-elements
% in Ht. These rows are different for every column of P (every row of Ht). The
% number of rows is, however, always equal to Ne+p
rr = (0:Ne-1) - floor((Ne-2)/2);
I = eye(Ne);
for j = 1:Ne,
   rows = midpoint(j)+rr;
   nn = find(rows<1);
   pp = find(rows>N);
   if periodic,
      rows(nn) = rows(nn)+N;
      rows(pp) = rows(pp)-N;
   else
      if ~isempty(nn),
         rows = (1:Ne);
      elseif ~isempty(pp),
         rows = (N-Ne+1:N);
      end
   end
   Prows = P(rows,1:Ne);
   b = I(j,1:Ne);
   % pseudo-inverse: system is often underdetermined
   % This assumes that a solution with zeros outside the selected rows exists. 
   % It is possible that the system is undertermined and that adding nonzeros 
   % outside these rows is the only way to find a solution for the original 
   % problem. We however assume that a solution exists within these rows, and
   % that the underdetermined system can be solved exactly with the 
   % pseudo-inverse. This is certainly the case if the bandwidth of the system
   % is smaller than the size of the vector "rows": going outside these rows 
   % would then add only zeros
   Htj = sparse((pinv(full(Prows'))*b')');
   Ht(j,rows) = Htj;

   % now we find one row of the homogeneous solution
   % We solve Prows0'*Htj0 = 0.
   % Put A = Prows0' and X = Htj0, and suppose X(i) = 1; i.e., i = midpoint(j)
   % then AX = 0 <=> A(:,not_i)*X(not_i) = -A(:,i)*X(i) = - A(:,i)
   rows0 = setdiff(rows,midpoint(j));
   Prows0 = P(rows0,1:Ne);
   b = -P(midpoint(j),1:Ne);
   Htj0 = sparse((pinv(full(Prows0'))*b')');
   Ht0(j,rows0) = Htj0;
   Ht0(j,midpoint(j)) = 1;
end

% now find Ht+A*Ht0 so that vanishing moment conditions are met
% Solve a system BXC = R with X bandlimited matrix, B = M'*P; C = Ht0;
if p > 0,
   I = eye(N);
   B = M'*P;
   C = Ht0;
   R = M'*(I-P*Ht);
   A = squaresparsesolveAXBisC(B,C,R);
   Ht = Ht + A*Ht0;
end
