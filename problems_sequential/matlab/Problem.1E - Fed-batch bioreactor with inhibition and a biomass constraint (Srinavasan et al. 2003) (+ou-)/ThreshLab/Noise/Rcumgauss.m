function Rcumgauss = Rcumgauss(x,s);

% Rcumgauss -- survival function --- complementary cumulative distribution 
%                                    function for zero-mean normal RV
%   Usage
%     f = Rcumgauss(x,s);
%   Inputs
%     x     value in which 1-cdf is evaluated
%     s     st.dev. of normal distribution (default 1)
%   Outputs
%     f     1-cdf
%   See also
%     help gauss (probability density function)
%     help cumgauss (cumulative distribution function)


if (nargin == 1)
   s = 1;
end

Rcumgauss = 0.5 * erfc(x./(sqrt(2)*s));

if s<eps,
   Rcumgauss(abs(x)<eps) = 0;
end


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
