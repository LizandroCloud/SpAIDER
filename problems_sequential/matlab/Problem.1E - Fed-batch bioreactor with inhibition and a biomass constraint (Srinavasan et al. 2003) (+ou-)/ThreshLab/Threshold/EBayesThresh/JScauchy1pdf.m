
function f = JScauchy1pdf(u);

% JScauchy1pdf -- probablity density function of Johnstone's and Silverman's
%                 Caucly-like prior
%   Usage
%     f = JScauchy1pdf(u);
%   Inputs
%     u     value in which pdf is evaluated
%   Outputs
%     f     pdf
%   Reference
%     I. Johnstone and B. Silverman, Needles and Straw in Haystacks, Annals of
%                                    Statistics, 2004.
%   See also
%     help JScauchypdf
%     help JScauchy1cdf
%     help JScauchyplusnormalpdf

u = abs(u);
freal = 1/sqrt(2*pi)-u.*exp(u.^2/2).*Rcumgauss(u,1);
freal = max(freal,0); % to eliminate minus infinities (due to instabilities)
u = max(u,1e-6); % to avoid division by 0
fapprox = 1./(sqrt(2*pi)*u.^2);
f = freal.*(u<37.5) + fapprox.*(u>=37.5);

