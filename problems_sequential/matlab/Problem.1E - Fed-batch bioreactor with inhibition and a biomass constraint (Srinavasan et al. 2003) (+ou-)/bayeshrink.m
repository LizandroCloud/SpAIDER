function T = bayeshrink(CD,CV)
[m,n] = size(CD);
sigma  = (median(median(abs(CD))))/0.6745; %Standard deviation
var_sig = std_sig^2;            % Variance
mean_sig = mean(mean(CH));      %mean value
sigmax = var_sig-(mean_sig^2);
sigmas = sqrt(max((sigmax-(sigma^2)),0));
%Calculation of threshold
if sigmas ~=0
    T = ((sigma)^2)/(sigmas);
elseif sigmas == 0
    T = max(max(abs(CV)));
end