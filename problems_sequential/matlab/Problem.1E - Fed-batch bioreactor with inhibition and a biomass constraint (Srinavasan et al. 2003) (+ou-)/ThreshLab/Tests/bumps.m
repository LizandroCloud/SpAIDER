
function f = bumps(x);

% ThreshLab/Tests/bumps -- creates bumps testfunction
%  Usage
%    f = bumps(x);
%    f = bumps(N);
%  Inputs
%    x      abscis value (possibly (non-)equidistant)
%    N      number of equidistant abscis values (1:N)/N
%  Outputs
%    f      function values of bumps testfunction in x or (1:N)/N
%  Description
%    evaluates bumps testfunction in t = (x - x(1))/(x(N)-x(1))
%  See also
%    help blocks
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
hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
f = zeros(size(t));
for j =1:length(pos)
   f = f + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
end

