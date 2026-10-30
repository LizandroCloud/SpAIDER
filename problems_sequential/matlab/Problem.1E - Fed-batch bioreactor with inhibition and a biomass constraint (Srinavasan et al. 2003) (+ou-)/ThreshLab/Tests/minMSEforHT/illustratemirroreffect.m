
% This routine illustrates the mirror effect in Cp and AIC when using in model
% selection in the presence of sparsity
%
% This routine plots results of SUREhard, MSEhard. See also GCVhard

n = 2000;
p = 0.05;
stdev = 1;

if exist('equidistantn1') ~= 1, equidistantn1 = 0; end
% The binary value equidistantn1 decides whether the plots are made for
% equidistant n1 (number of selected components) or equidistant threshold
% values.
m = 1000; % number of values n0 for which AIC/Cp is evaluated

if exist('realrandom') ~= 1, realrandom = 0; end

if ~realrandom,
   randn('state',2000);
   rand('state',2000);
end



omega = randn(1,n)*stdev;

a = 1/5;
r = rand(1,n);
v = randlaplace(1,n,a).*(r<p);
[vsort i] = sort(abs(v)); vsort = v(i); v = vsort(length(v):-1:1);


w = v+omega;

W = w;


X = W.^2; [X s] = sort(X); W = W(s);
alfa = (1:n)/n;
lambda = invcumgauss(1/2+alfa/2,stdev);

if equidistantn1,
   s = floor(n/m);
   nn0 = (1:s:n); % vector of values n0 = n-n1; n1 is number of selected
                  % components
   m = length(nn0); 
   % Rounding effects may cause length(nn0) to be slightly different from the
   % original m
   t = abs(W(nn0)); % empirical threshold corresponding to values of n0
else
   tmax = max(abs(W));
   t = (1:m)/m*tmax;
   nn0 = zeros(1,m);
   for k = 1:m, nn0(k) = sum(abs(w)<t(k)); end
end
nn1 = n - nn0;

gcvsoft = GCV1soft(w,t);
gcvhard = GCV1softhard(w,t);
msesoft = MSE1(v,w,t,1);
msehard = MSE1(v,w,t,0);
STrisk = softrisk(t,v,stdev);
HTrisk = hardrisk(t,v,stdev);

suresoft = SURE1soft(w,t,stdev);
surehard = SURE1hard(w,t,stdev);

% surehardsoft is min(Cp; n1 components) as function of n1
surehardsoft = zeros(size(surehard));
surehardsoft1 = surehardsoft;
for k = 1:m
   wthr = HT(w,t(k));
   N1 = nn1(k);
   surehardsoft(k) = norm(wthr-w)^2/n +2*N1/n*stdev^2 - stdev^2;
   surehardsoft1(k) = norm(wthr-w)^2/n;
end
correction = 4*stdev^2*t.*gauss(t,stdev);

% Next is linear selection on oracle-ordered data
Cp = zeros(1,m); SPE = zeros(1,m); % SPE = squared prediction error
for k = 1:m
   N1 = nn1(k);
   wlin = zeros(size(w));
   wlin(1:N1) = w(1:N1);
   Cp(k) = norm(wlin-w)^2/n +2*N1/n*stdev^2 - stdev^2;
   SPE(k) = norm(wlin-v)^2/n;
end

figure(1)
plot(nn1,SPE,'m','linewidth',2)
hold on
plot(nn1,Cp,'b','linewidth',2)
plot(nn1,surehardsoft,'k','linewidth',2)
plot(nn1,msehard,'g','linewidth',2)
plot(nn1,surehardsoft+correction,'r--','linewidth',2)
legend('SPE oracle','Cp','surehardsoft','surehard','msehard')
hold off
xlabel('n_1: number of components in model','fontsize',14,'fontweight','b')

% Now the same for AIC vs R (Kulback-Leibler)

AIClin = zeros(1,m); AICthr0 = zeros(1,m);
Rlin = zeros(1,m); Rthr = zeros(1,m);
sigma2hatthr = zeros(1,m);
sigma2hatlin = zeros(1,m);
for k = 1:m
   N1 = nn1(k);
   wlin = zeros(size(w));
   wlin(1:N1) = w(1:N1);
   sigma2hatlin(k) = norm(w-wlin)^2/(n-N1);
   AIClin(k) = -(n-N1)/n - log(2*pi*sigma2hatlin(k)) - 2*N1/n;
   Rlin(k) = -norm(wlin-v)^2/sigma2hatlin(k)/n - log(2*pi*sigma2hatlin(k)) ...
             -stdev^2/sigma2hatlin(k);

   wthr = HT(w,t(k));
   sigma2hatthr(k) = norm(w-wthr)^2/(n-N1);
   AICthr0(k) = -(n-N1)/n - log(2*pi*sigma2hatthr(k)) - 2*N1/n;
   Rthr(k) = -norm(wthr-v)^2/sigma2hatthr(k)/n - log(2*pi*sigma2hatthr(k)) ...
             -stdev^2/sigma2hatthr(k);
end

stdevestGCV = sqrt(varestGCV(w));
n0tildestGCV = n*(cumgauss(t,stdevestGCV)-cumgauss(-t,stdevestGCV));
AICthrcorrGCV = -(n-nn1)/n - log(2*pi*sigma2hatthr) ...
              -2*n0tildestGCV/n.*(nn1/n+2*t.*gauss(t,stdevestGCV))./...
                                 (n0tildestGCV/n-2*t.*gauss(t,stdevestGCV));

stdevestMAD = MAD1(w)/invcumgauss(3/4);
n0tildestMAD = n*(cumgauss(t,stdevestMAD)-cumgauss(-t,stdevestMAD));
AICthrcorrMAD = -(n-nn1)/n - log(2*pi*sigma2hatthr) ...
              -2*n0tildestMAD/n.*(nn1/n+2*t.*gauss(t,stdevestMAD))./...
                                 (n0tildestMAD/n-2*t.*gauss(t,stdevestMAD));


figure(2)
plot(nn1,Rlin,'m','linewidth',2)
hold on
plot(nn1,AIClin,'b','linewidth',2)
plot(nn1,AICthr0,'k','linewidth',2)
plot(nn1,Rthr,'g','linewidth',2)
plot(nn1,AICthrcorrGCV,'r--','linewidth',2)
plot(nn1,AICthrcorrMAD,'m--')
legend('R oracle','AIC oracle','AICthr','Rthr','AIC with mirror corr.')
hold off
xlabel('n_1: number of components in model','fontsize',14,'fontweight','b')
axis([0,n,-10,1])


