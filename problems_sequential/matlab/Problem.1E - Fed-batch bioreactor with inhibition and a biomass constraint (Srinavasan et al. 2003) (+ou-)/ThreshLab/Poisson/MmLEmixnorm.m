
function p = MmLEmixnorm(w,sig0,sig1);

% usage: p = MmLEmixnorm(w,sig0,sig1);
% maximum marginal likelihood estimation of parameter p in mixed normal model:
% p = P(coefficient has normal density with st.dev = sig0)
% 1-p = P(coefficient has normal density with st.dev = sig1)
% w: observed values
% sig0 and sig1 are standard deviations

eps = 0.001; %i.e. 10 function evaluations
p0 = 0.0;
p1 = 1.0;
interval = [p0 p1];
L0grad = dlnLmixednormaldp(p0,w,sig0,sig1);
L1grad = dlnLmixednormaldp(p1,w,sig0,sig1);
if (L0grad*L1grad) > 0
   L0 = lnLmixednormal(p0,w,sig0,sig1);
   L1 = lnLmixednormal(p1,w,sig0,sig1);
   p = double(L1 > L0);
else
  while (p1-p0) > eps
    p = (p0+p1)/2;
    Lgrad = dlnLmixednormaldp(p,w,sig0,sig1);
    if ((Lgrad*L0grad) > 0)
        p0 = p;
        L0grad = Lgrad;
    else
        p1 = p;
        L1grad = Lgrad;
    end
  end
p = (p0+p1)/2;
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
