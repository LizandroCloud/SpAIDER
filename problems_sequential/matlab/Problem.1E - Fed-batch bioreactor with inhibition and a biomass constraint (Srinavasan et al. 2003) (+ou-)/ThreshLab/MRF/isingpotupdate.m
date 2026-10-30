
function Uc = isingpotupdate(X0,X1,tau)

% Uc = isingpotupdate(X0,X1,tau)
% X0 and X1 are matrices with elements 1 and/or -1
% Uc contains the potential values when entries of X1 are surrounded by
% corresponding neighboring pixels in X0
% This can be used to examine in which pixels a switch in label values
% (1 <-> -1) would be advantageous from the energetic point of view.

[N M] = size(X0);
positions = [-1 -1
             -1  0
             -1  1
              0 -1
              0  0
              0  1
              1 -1
              1  0
              1  1];
Uc = zeros(size(X0(2:N-1,2:M-1)));
Ns = 0;
for s = 1:9
   sx = positions(s,1);
   sy = positions(s,2);
   if (s == 5)
     Uc = Uc + 0;
   else
     Uc = Uc - X0(2+sx:N-1+sx,2+sy:M-1+sy).*X1(2:N-1,2:M-1);
   end
end
% H = tau*sum(sum(Uc));
Uc = tau*Uc;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
