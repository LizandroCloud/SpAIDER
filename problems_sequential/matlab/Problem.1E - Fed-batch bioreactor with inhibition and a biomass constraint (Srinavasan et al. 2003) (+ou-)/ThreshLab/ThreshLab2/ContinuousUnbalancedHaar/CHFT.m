
function z = CHFT(y,M,m);

% CHFT - ThreshLab2/ContinuousUnbalancedHaar/ -- Continuous Haar Fisz Transform
%   Usage
%     z = CHT(y,M,m);
%   Inputs
%     y   observed values of function to be analysed
%     M   number of dyadic levels/scales (depth of transform)
%     m   number of intermediate scales (see below; default is 0)
%   Outputs
%     z   continuous wavelet transform coefficients
%   Description
%     evaluates the Continuous Haar Fisz Transform in a finite number of scales,
%     specified by parameters M and m. M is the number of dyadic scales, as in
%     the Discrete Wavelet Transform (decimated or undecimated).
%     m is number of additional, intermediate scales; m is ceiled to the next
%     value of 2^(r+1)-r-2 for integer r. This is to make an equidistant
%     distribution of the additional scales in between the coarsest dyadic
%     scales.
%     If m = 0, then this function coincides with the undecimated Haar wavelet
%     transform, see  RWT_BIO
%   See also
%      help RWT_BIO
%      help FWT_BIO
%      help CHT
%      help FFHFW

[w s] = CHT(y,M,m);
s = s + double(s==0);
z = w./sqrt(s);

