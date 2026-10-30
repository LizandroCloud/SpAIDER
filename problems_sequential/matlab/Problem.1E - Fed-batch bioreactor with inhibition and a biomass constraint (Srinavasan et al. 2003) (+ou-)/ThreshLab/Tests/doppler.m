function f = doppler(x);
                                                                                
% doppler -- creates doppler testfunction
%  Usage
%    f = doppler(x);
%    f = doppler(N);
%  Inputs
%    x      abscis value (possibly (non-)equidistant)
%    N      number of equidistant abscis values (1:N)/N
%  Outputs
%    f      function values of doppler testfunction in x or (1:N)/N
%  Description
%    evaluates doppler testfunction in t = (x - x(1))/(x(N)-x(1))
%  See also
%    help blocks
%    help bumps
%    help heavisine
%    help skyline

N = length(x);
if N == 1,
   N = x;
   x = (1:N);
end
t = (x - x(1))/(x(N) - x(1));
f = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05));
