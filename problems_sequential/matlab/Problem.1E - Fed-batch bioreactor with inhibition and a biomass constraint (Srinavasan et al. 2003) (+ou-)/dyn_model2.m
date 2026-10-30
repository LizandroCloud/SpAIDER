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


dx(1) = {'( mim*x(2) )'};
dx(2) = {'- ( mim*x(2) )'};


end