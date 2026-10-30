function [ dx ] = sens( t, x, w, ks )
% State system
dx=[ x(2) + w(ks);
- w(ks);
x(1)^2 / 2. ];
% Append sensitivity system
for is = 1:length(w)
    if is == ks
uws = 1.;
 else
 uws = 0.;
 end
 dx = [ dx;
 x(3*is+2) + uws;
 - uws;
 x(1) * x(3*is+1) ];
 end
 end