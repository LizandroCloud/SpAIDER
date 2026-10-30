function gauss = gauss(x,s)

% gauss -- probability density function (pdf) for zero-mean normal RV
%   Usage
%     f = gauss(x,s);
%   Inputs
%     x     value in which pdf is evaluated
%     s     st.dev. of normal distribution (default 1)
%   Outputs
%     f     pdf
%   See also
%     help cumgauss (cumulative distribution function)
%     help Rcumgauss (survival function)


if (nargin == 1)
   s = 1;
end

if s>eps,
   gauss = exp(-x.*x./(2*s.*s))./(sqrt(2*pi)*s);
else
   gauss = zeros(size(x));
   gauss(abs(x)<eps) = Inf;
end


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
