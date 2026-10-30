function [value, isterminal, direction]=eventfcn(t,y)
global cAf;
value = y(3)-cAf;
isterminal = 1;
direction = 0;