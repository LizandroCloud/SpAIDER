function [h, ht] = Haar();

%   Usage
%     h = Haar;
%     [h ht] = Haar();
%   Inputs
%     none
%   Outputs
%     h    Haar filter coefficient
%     ht   dual Haar filter coefficient = h (orthogonal transform)
%   Description
%     creates Haar filter
%   Note
%     This file called by makeWTfilter
%   See also
%     help makeWTfilter

h = [1 1]/sqrt(2);
ht = h;
