
function f = heavisine(x);

% ThreshLab/Tests/heavisine -- creates heavisine testfunction
%  Usage
%    f = heavisine(x);
%    f = heavisine(N);
%  Inputs
%    x      abscis value (possibly (non-)equidistant)
%    N      number of equidistant abscis values (1:N)/N
%  Outputs
%    f      function values of heavisine testfunction in x or (1:N)/N
%  Description
%    evaluates heavisine testfunction in t = (x - x(1))/(x(N)-x(1))
%  See also
%    help blocks
%    help bumps
%    help doppler
%    help skyline

N = length(x);
if N == 1,
   N = x;
   x = (1:N);
end
t = (x - x(1))/(x(N) - x(1));
f = 4.*sin(4*pi*t);
f = f + 2 - 2*double(t>.3) - 2*double(t<.72);
   % note that this definition is a bit different from WaveLab, where the
   % discontinuities were defined with the sign-function. This sign function
   % has values -1,0, and 1, which may lead to three different values before,
   % in and after a discontinuity.
   % f = f - sign(t - .3) - sign(.72 - t);
