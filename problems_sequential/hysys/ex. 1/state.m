function [ dx ] = state( t, x, w, ks )
% State system
dx=[ x(2) + w(ks);
- w(ks);
x(1)^2 / 2. ];
end