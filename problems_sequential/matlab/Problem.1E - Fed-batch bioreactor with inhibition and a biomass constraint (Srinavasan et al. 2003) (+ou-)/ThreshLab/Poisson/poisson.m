
function y = poisson(lambda,realrandom);

% usage: y = poisson(lambda,realrandom);
% or: y = poisson(x); 
% default for realrandom is 0.
% generates poisson-noise with intensities lambda.

if nargin < 2 | ~realrandom,
   rand('seed',931316785);
end
y = zeros(size(lambda));
for i = 1:size(lambda,1),
   for j = 1:size(lambda,2),
     y(i,j) = 0;
     sumarrivaltime = -log(rand);
     while sumarrivaltime < lambda(i,j)
       sumarrivaltime = sumarrivaltime - log (rand);
       y(i,j) = y(i,j) + 1;
     end
   end
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
