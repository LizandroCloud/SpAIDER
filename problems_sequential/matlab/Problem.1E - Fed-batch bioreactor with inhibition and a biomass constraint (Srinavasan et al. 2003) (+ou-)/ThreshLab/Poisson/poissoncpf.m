
function p = poissoncpf(x,lambda,U);

% usage: p = poissoncpf(x,lambda);
% Poisson Cumulative Probability Function

if nargin<3,
   U = 100; % upperbound for exact computation, above U, the normal
            % approximation is used.
end
if any(size(lambda) ~= size(x)),
   if length(lambda) == 1,
      lambda = lambda*ones(size(x));
   elseif length(x) == 1,
      x = x*ones(size(lambda));
   else
      error('input vector sizes don''t match')
   end
end
[N, M] = size(x);
x = floor(x);
p = zeros(size(x));
for i = 1:N, for j = 1:M,
    if x(i,j) < U,
       p(i,j) = 0;
       for k = 0:x(i,j),
          p(i,j) = p(i,j) + lambda(i,j)^k/factorial(k);
       end;
       p(i,j) = exp(-lambda(i,j))*p(i,j);
    else
       z = (x(i,j)+0.5-lambda(i,j))/sqrt(lambda(i,j));
       p(i,j) = cumgauss(z);
    end
end, end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
