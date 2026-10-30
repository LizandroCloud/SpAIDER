
% Simulation study routine for GCV based robust variance estimation.
% This simulation study compares varestGCV with median absolute deviation (MAD)
% Starting the routine will display information on parameters
%
%
% See also 
%    help varestGCV

disp('This routine verifies the robustness of a variance estimator using')
disp('residuals from GCV thresholding')
disp('It uses the following variables (with default values)')
disp('Re-run after setting your own values if necessary')
disp('------------------------------------------------------------------')
disp('n        number of observations per simulation; default =  100')
disp('nsimul   number of simulations;                            200')
disp('m        length of vector pp                                50')
disp('pmax     maximum value of vector pp                        1.0')
disp('pp = (1:m)/m*pmax  proportions of outliers (number of nonzeros in V)')
disp('stdev    standard deviation of noise U in observational model W = V+U')
disp('                                                 default =   1')
disp('priormodel model for outliers in V: ''normal'',''Cauchy'',''Bernoulli''')
disp('           or ''Laplace'' (=double exponential) default = ''Bernoulli''')
disp('laplacepdfparameter  parameter in Laplace prior model')

disp('Press <Enter> to continue or <Crtl>C to stop and set other values')
disp('type >> pause off     before re-running for uninterrupted execution')
pause


if exist('n') ~= 1, n =  100; end
if exist('nsimul') ~= 1, nsimul = 200; end
if exist('m') ~= 1, m = 50; end
if exist('pmax') ~= 1, pmax = 1.0; end
if exist('stdev') ~= 1, stdev = 1; end
if exist('laplacepdfparameter') ~= 1, laplacepdfparameter = 1/5; end
if exist('priormodel') ~= 1, priormodel = 'Bernoulli'; end
rand('state',0)
randn('state',0)

shatMAD = zeros(nsimul,m);
shatGCV = zeros(nsimul,m);
pp = (1:m)/m*pmax;

for k = 1:m,
   k
   p = pp(k);
   for simul = 1:nsimul,
      r = rand(1,n);
      label = find(r<p);

      v = zeros(1,n);
      
      switch priormodel(1:2)
      case {'la','La','de','De','DE'}, % Laplace, double-exponential
         a = laplacepdfparameter;
         v(label) = randlaplace(1,size(label),a);
      case {'ca','Ca'}, % Cauchy
         v(label) = randcauchy(size(label));
      case {'no'}, % normal
         v(label) = randn(size(label))*stdev*6;
      case {'bi','be','Be'} % binary, Bernoulli
         v(label) = sign(rand(size(label))-0.5)*stdev*10;
      otherwise
         v(label) = 0;
      end
      w = v+randn(size(v))*stdev;
      shatMAD(simul,k) = MAD1(w)/invcumgauss(3/4);
      shatGCV(simul,k) = sqrt(varestGCV(w,1));
   end
end

shatGCVmean = mean(shatGCV);
shatMADmean = mean(shatMAD);

shatGCVmedian = median(shatGCV);
shatMADmedian = median(shatMAD);

shatGCVmax = max(shatGCV);
shatMADmax = max(shatMAD);

shatGCVmin = min(shatGCV);
shatMADmin = min(shatMAD);




figure(1)
plot(pp,shatMADmean,'m--','linewidth',2)
hold on
plot(pp,shatGCVmean,'k','linewidth',2)
hold off
legend('MAD','GCV','location','NW')
title('Mean values of stdev. estimators')
xlabel('proportion of outliers')
ylabel('estimator of st.dev.')

figure(2)
plot(pp,shatMADmedian,'m--','linewidth',2)
hold on
plot(pp,shatGCVmedian,'k','linewidth',2)
hold off
legend('MAD','GCV','location','NW')
title('Median values of stdev. estimators')
xlabel('proportion of outliers')
ylabel('estimator of st.dev.')


figure(3)
plot(pp,stdev./shatMADmin,'m--','linewidth',2)
hold on
plot(pp,stdev./shatGCVmin,'k','linewidth',2)
hold off
legend('MAD','GCV','location','NW')
title('inverse min. values of stdev. estimators')
xlabel('proportion of outliers')
ylabel('estimator of st.dev.')

figure(4)
plot(pp,shatMADmax./stdev,'m--','linewidth',2)
hold on
plot(pp,shatGCVmax./stdev,'k','linewidth',2)
hold off
legend('MAD','GCV','location','NW')
title('Max values of stdev. estimators')
xlabel('proportion of outliers')
ylabel('estimator of st.dev.')

