
function f = makesignal(x,number)

% usage: f = makesignal(x,number)
% or   : f = makesignal(N,number)
% or   : f = makesignal(N,name)
% or   : f = makesignal(x,name)
%
% input:
% x: grid
% N: number of equidistant grid points on [0,1]
% number: number of testfunction
% name  : name of testfunction
%          'heavisine','heavysine','block(s)','bump(s)','doppler','skyline'


if ischar(number) | isinf(number(1)),
   name = number;
   number = Inf;
end

N = length(x);
if N == 1,
   N = x;
   x = (1:N)/N;
end

if (number == 0)
   f = zeros(size(x));
elseif (number == 1)
   f = ones(size(x));
elseif (number == 2)
   f = sin(3*pi*x/4);
elseif (number == 3)
   f = skyline(x);
elseif (number == 4)
   % Donoho's heavy sine 
   f = 4.*sin(4*pi.*x);
   f = f - sign(x - .3) - sign(.72 - x);
elseif (number == 401)
   % Donoho's heavy sine on [0,1]
   t = (x - x(1))/(x(N) - x(1));
   f = 4.*sin(4*pi.*t);
   f = f - sign(t - .3) - sign(.72 - t);
elseif (number == 5)
   f = sin(4*pi.*x).*exp(-x);
elseif (number == 51)
   f = sin(4*pi.*x).*exp(-x);
   f = f + (x>0.3) - (x>0.7)*0.3;
elseif (number == 6)
   f = sin(12.*x).*exp(-x);
elseif (number == 61)
   f = sin(4*pi.*x).*exp(-x);
   f = f + (x>0.3) - (x>0.7)*0.3;
else
   f = ones(size(x));
end

if isinf(number)
switch name(2:4),
   case 'eav',
      f = heavisine(x);
   case 'loc',
      f = blocks(x);
   case 'ump',
      f = bumps(x);
   case 'opp',
      f = doppler(x);
   otherwise
      f = skyline(x);
end, end
