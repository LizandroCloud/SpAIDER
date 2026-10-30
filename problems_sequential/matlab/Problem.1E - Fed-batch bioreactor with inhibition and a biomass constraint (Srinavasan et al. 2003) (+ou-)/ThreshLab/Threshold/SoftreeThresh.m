
% usage:
% wt = SoftreeThresh(w,t);
% or:
% [wt, label] = SoftreeThresh(w,t);

function [wt, label] = SoftreeThresh(w,t);

wt = w;
label = ones(size(w));
n = length(w);
J = 0;
while (2^J < n)
   J = J+1;
end;
L = log2(2^J-n);
M = J-L;
for j = L:J-2
for i = 2^j+1-2^L:2^(j+1)-2^L
  if label(i)
     if (abs(w(i)) < t)
        label(i) = 0;
     else
        wt(i) = sign(w(i))*(abs(w(i)) - t);
     end
  end
  if (label(i) == 0)
     wt(i) = 0;
     l = 2*i + 2^L - 1;
     label(l:l+1) = 0;
  end
end   
end
for i = 2^(J-1)+1-2^L:2^J-2^L
  if label(i)
     if (abs(w(i)) < t)
        label(i) = 0;
     else
        wt(i) = sign(w(i))*(abs(w(i)) - t);
     end
  end
  if (label(i) == 0)
     wt(i) = 0;
  end
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
