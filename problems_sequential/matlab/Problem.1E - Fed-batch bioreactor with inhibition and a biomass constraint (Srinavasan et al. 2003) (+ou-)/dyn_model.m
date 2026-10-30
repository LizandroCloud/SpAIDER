function [ dx par] = dyn_model

% ************************************************************************
%                              #Case study C#
% Problem:     Fed-batch bioreactor with inhibition and a biomass constraint  
% References:  Srinivasan et al. (2003a) 
%              Visser, Srinivasan, Palanki, & Bonvin, (2000)
%
% ************************************************************************

par(1)={'mim = 0.53'}; % l/(h)
par(2)={'km = 1.2'};   % (g/L)
par(3)={'ki = 22'}; % (g/L)
par(4)={'Yx = 0.4'};
par(5)={'Yp = 1'};
par(6)={'v = 0.5'}; %L/h;
par(7)={'Sin = 20'}; % g/l
par(8)={'xo = 1'}; %g/L;
par(9)={'So = 0'}; %g/L;
par(10)={'Po = 0'}; %g/L;
par(11)={'Vo = 2'}; %L

dx(1) = {'( mim*x(2) / (km + x(2)+ ( (x(2)^2)/ki )) )*x(1) - u*x(1)/x(4)'};
dx(2) = {'- ( mim*x(2) / (km + x(2)+ ( (x(2)^2)/ki )) )*x(1)/Yx - v*x(1)/Yp + u*(Sin-x(2))/x(4)'};
dx(3) = {'v*x(1) - u*x(3)/x(4)'};
dx(4) = {'u'};

end