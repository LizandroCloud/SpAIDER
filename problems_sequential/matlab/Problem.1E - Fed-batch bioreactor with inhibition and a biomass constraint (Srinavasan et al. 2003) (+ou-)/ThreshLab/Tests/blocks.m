
function f = blocks(x);

% ThreshLab/Tests/blocks -- creates blocks testfunction
%  Usage
%    f = blocks(x);
%    f = blocks(N);
%  Inputs
%    x      abscis value (possibly (non-)equidistant)
%    N      number of equidistant abscis values (1:N)/N
%  Outputs
%    f      function values of blocks testfunction in x or (1:N)/N
%  Description
%    evaluates blocks testfunction in t = (x - x(1))/(x(N)-x(1))
%  See also
%    help bumps
%    help doppler
%    help heavisine
%    help skyline

N = length(x);
if N == 1,
   N = x;
   x = (1:N);
end
t = (x - x(1))/(x(N) - x(1));
pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
f = zeros(size(t));
for j=1:length(pos)
   f = f + double(t >= pos(j)).*hgt(j) ;
   % note that this definition is a bit different from WaveLab, where the
   % discontinuities were defined with the sign-function. This sign function
   % has values -1,0, and 1, which may lead to three different values before,
   % in and after a discontinuity.
   % f = f + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
end
