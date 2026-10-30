
function tau = MPLEising(X);

% usage: tau = MPLEising(X);
% Maximum Pseudo Likelihood Estimation for rigity parameter in Ising model

maxtau = 1;
eps = 0.0001;
last = maxtau;
first = 0;

F0 = 1;
F1 = 1;
F2 = 2;
steps = 1;

while ((last - first)/F2 > (eps* last))
   F0 = F1;
   F1 = F2;
   F2 = F0 + F1;
   steps = steps + 1;
end
disp('MPLE: ')
steps
v = first + (F1/F2) * (last-first);
fv = PLising(X,v);
u = first + (F0/F2) * (last-first);
fu = PLising(X,u);
for step = 1:steps
   if (fu < fv)
      first = u;
      u = v;
      fu = fv;
      v = last - u + first;
      fv = PLising(X,v);
   else
      last = v;
      v = u;
      fv = fu;
      u = first + last - v;
      fu = PLising(X,u);
   end
end
tau = (u+v)/2;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
