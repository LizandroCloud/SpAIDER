function x = logfactorial(k);

% logfactorial - ThreshLab/General/ -- compute log(k!) for integer k
%   Usage
%     n = logfactorial(k);
%   Inputs
%     k integer (or vector of integers)
%   Outputs
%     x = log(k!)
%   Description
%     Uses Stirling's approximation k! ~ sqrt(2*pi*k) * k^k/e^k for large
%     values of k

x = zeros(size(k));
smallk = find(k<=40);
largek = find(k>40);

k1 = k(smallk);
k1 = max(k1,1);
k1fac = k1;
while max(k1) > 1, k1 = max(k1-1,1); k1fac = k1fac.*k1; end
x(smallk) = log(k1fac);

k = k(largek);
x(largek) = k.*log(k)-k+log(k)/2+log(2*pi)/2+1./(12*k)-1./(360*k.^3);
