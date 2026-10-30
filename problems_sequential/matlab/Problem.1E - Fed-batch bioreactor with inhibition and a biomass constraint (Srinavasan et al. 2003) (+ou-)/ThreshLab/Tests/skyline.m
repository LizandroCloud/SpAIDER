
function f = skyline(x);

% ThreshLab/Tests/skyline -- creates skyline testfunction
%  Usage
%    f = skyline(x);
%    f = skyline(N);
%  Inputs
%    x      abscis value (possibly (non-)equidistant)
%    N      number of equidistant abscis values (1:N)/N
%  Outputs
%    f      function values of skyline testfunction in x or (1:N)/N
%  Description
%    evaluates skyline testfunction in x/max(x)
%  See also
%    help blocks
%    help bumps
%    help doppler
%    help heavisine

N = length(x);
if N == 1,
   N = x;
   x = (1:N);
end
% t = (x - x(1))/(x(N) - x(1));
% t = 2048*t;

t = x*2048/x(N);
f = (t <= 128)*2;
f = f+(t>128).*(t<=160)*5;
f = f+(t>160).*(t<=320);
f = f+(t>336).*(t<=512)*3;
f = f+(t>512).*(t<=544)*7;
f = f+(t>=800).*(t<=960).*(t-800)/160;
f = f+(t>960).*(t<=1024);
f = f+(t>1024).*(t<=1200)*3;
f = f+(t>1200).*(t<=1360)*2;
f = f+(t>1360).*(t<=1600).*(240-t+1360)/240;
f = f+(t>1615).*(t<1632).*10;
f = f+(t>=1632).*(t<=1760).*(t-1632)/128;
f = f+(t>1760);
