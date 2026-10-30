
function [Y, K, x, stdev, Z] = ...
    setupsimulationYisKxplussigmaZ(typeofK,m,varargin);

% setupsimulationYisKxplussigmaZ -- set up test with Y = Kx + sigma*Z
%  Usage
%    [Y K x stdev Z] =  setupsimulationYisKxplussigmaZ(typeofK,m,n,p,...
%                       snratio,realrandom,
%  Inputs
%    typeofK  possible models for design matrix. Currently available:
%             'random' random matrix
%             'randomtoeplitz' random selection of n rows from m x m 
%                              stationary filter matrix
%             'randombandtoeplitz' random selection of n rows from m x m 
%                                  bandlimited stationary filter matrix
%             'randomfromFourier' random selection of n rows from m x m Fourier
%                                 transformation matrix
%             'SVDgivenSV' U'*Sigma*V, with U and V random orthogonal matrices
%                          and Sigma = diag(sigma), where sigma is given as
%                          optional parameter in 'singularvalues'
%             'SVDrandomSV' U'*Sigma*V, with U and V random orthogonal matrices
%                          and Sigma = diag(sigma), where sigma r random
%                          singular values, where r given as optional parameter
%                          in 'n_singval'
%             'Vandermonde' Vandermonde matrix for given or random covariate
%                           values (declared as ,'vdmcovar',vdmcovar)
%    m      model size of x
%    n      number of observations; default n = m
%    p      sparsity parameter in zero-inflated Laplacian model of x:
%           p = P(X =/= 0)      default p = 0.05
%    snratio signal-to-noise ratio; SNR = 20 log(norm(K*x)/norm(noise))
%            default = 10;
%    realrandom  if 0, then seed for random number generators is fixed
%                if 1, then seed is current value of seed
%                default is 0
%    'stdev',stdev use value stdev as standard deviation of nosie instead of
%                  signal-to-noise ratio
%    
%  Outputs
%    Y = K*x+stdev*Z
%    Y      observations, vector of size n x 1
%    K      design matrix, matrix of size n x m
%    x      unobserved model, vector of size m x 1, generated according to zero
%           inflated Laplacian model with a = 1/5 and p from input
%    stdev  standard deviation of noise
%    Z      standardized observational noise, IID, normal vector of size n x 1
%  Examples of usage
%
%    use
%       x = setupsimulationYisKxplussigmaZ('eye',n,n,p,Inf,realrandom);
%    to generate sparse vector only (no K, no errors)
%
%    if exist('typeofK') ~= 1, typeofK = 'random'; end
%    if exist('realrandom') ~= 1, realrandom = 0; end
%    if exist('dimK') ~= 1, dimK = 0; end
%       if ~dimK, m = 1000; n = 200; end % n = number of observations
%    clear dimK
%    if exist('snratio') ~= 1, snratio = 10; end
%    if exist('degreeofsparsity') ~= 1, degreeofsparsity = 0.05; end
%    if exist('bandwidth') ~= 1, bandwidth = 1; end
%    if exist('singvalues') ~= 1, singvalues = min(m,n); end
%    if exist('n_singval') ~= 1, n_singval = 1; end
%    if exist('standardize') ~= 1, standardize = 0; end
%    if exist('vdmcovar') ~= 1, vdmcovar = sort(rand(1,n)); end
%
%    p = degreeofsparsity;
%
%    [w K v stdev Z] = ...
%      setupsimulationYisKxplussigmaZ(typeofK,m,n,p,snratio,realrandom,...
%                                    'bandwidth',bandwidth,...
%                                    'singularvalues',singvalues,...
%                                    'n_singval',n_singval,...
%                                    'vdmcovar',vdmcovar);
%    
% 

if nargin < 2,
   error('typeofK and m are mandatory inputs')
end
nnargin = nargin;
ar = 3; 
while ar < nargin, 
   if ischar(varargin{ar-2}), 
      nnargin = ar-1; 
      ar = nargin;
   end
   ar = ar+1;
end
if nnargin < 3, n = m; else n = varargin{1}; end
if nnargin < 4, p = 0.05; else p = varargin{2}; end
if nnargin < 5, snratio = 10; else snratio = varargin{3}; end
if nnargin < 6, realrandom = 0; else realrandom = varargin{4}; end

nnargin = nnargin+1;
stdev0 = NaN;
vdmcovar = sort(rand(1,n));
r = min(m,n);
for ar = nnargin:nargin, arg = varargin{ar-2}; if ischar(arg)
  switch arg
  case {'bandwidth'}
    bandwidth = varargin{ar-1};
  case {'singularvalues'}
    sigma = varargin{ar-1};
  case {'n_singval'}
    r = varargin{ar-1};
  case {'stdev'}
    stdev0 = varargin{ar-1};
  case {'vdmcovar'}
    vdmcovar = varargin{ar-1};
  end
end, end
if typeofK(1:3) == 'SVD',
   if typeofK(4:6) == 'giv',
      r = length(sigma);
   else
      sigma = abs(rand(1,r));
   end
end
if length(typeofK)>6, if typeofK(2:7) == 'anderm',
   n = length(vdmcovar);
end, end

disp('-------------------------------------------------')
disp('Generating Y = K*x+stdev*Z with following values:')
disp(['   typeofK = ' typeofK])
disp(['         n = ' num2str(n)])
disp(['         m = ' num2str(m)])
disp(['         p = ' num2str(p)])
disp(['       snr = ' num2str(snratio)])
disp(['realrandom = ' num2str(realrandom)])
disp('-------------------------------------------------')
disp(' '), disp(' ')

if ~realrandom,
   randn('state',2000);
   rand('state',2000);
end

Keye = 0;
switch typeofK,
case {'randomtoeplitz'}
   h = randn(1,m);
   K = toeplitz(h);
   r = randperm(m);
   r = sort(r(1:n));
   K = K(r,1:m);
case {'randombandtoeplitz'}
   h = randn(1,m);
   h(bandwidth+1:m) = 0;
   K = toeplitz(h);
   r = randperm(m);
   r = sort(r(1:n));
   K = K(r,1:m);
case {'random'}
   K = randn(n,m);
case {'randomband','bandrandom'}
   dotranspose = 0;
   if n>m, dotranspose = 1; nm = n; n = m; m = nm; end
   K = sparse(diag(randn(1,m)));
   for d = 1:bandwidth, 
       K = K+sparse(diag(randn(1,m-d),d))+sparse(diag(randn(1,m-d),-d));
   end
   r = [0 floor((1:m)/m*n)]; r = find(r(2:m+1)~=r(1:m));
   K = K(r,1:m);
   if dotranspose == 1, K = K'; nm = n; n = m; m = nm; end
case {'eye'}
   n = m;
   K = eye(m);
   Keye = 1;
case {'randomfromFourier'}
   halfm = floor(m/2-1/2);
   restm = m-1-halfm*2;
   K = fft(eye(m));
   K = K(1:m,1:1+halfm+restm);
   K = [real(K(:,1:1+halfm+restm)) imag(K(:,2:1+halfm))];
   K = K';
   r = randperm(m);
   r = sort(r(1:n));
   K = K(r,1:m);
case {'SVDgivenSV','SVDrandomSV'}
   mn = min(m,n);
   r = min(mn,r); sigma = sigma(1:r);
   sigma0 = zeros(1,mn); sigma0(1:r) = sigma;
   Sigma0 = diag(sigma0);
   Sigma = zeros(n,m); Sigma(1:mn,1:mn) = Sigma0;
   U = randn(n,n); V = randn(m,m);
   U = gramschmidt(U); V = gramschmidt(V);
   K = U'*Sigma*V;
case {'Vandermonde'}
   K = zeros(n,m);
   for k=1:m,
      K(1:n,k) = column(vdmcovar.^(k-1));
   end
end
if ~Keye, % binary variable Keye saves computation time in case K = eye(n)
   K = (1/norm(full(K)))*K;
end

a = 1/5;
r = rand(m,1);
x = randlaplace(m,1,a).*(r<p);


if ~Keye, Kx = K*x; else, Kx = x; end
if isinf(snratio), Y = Kx; else, Y = makenoisy(Kx,snratio); end
U = Y-Kx;
stdev = sqrt(var2(U));
if ~isnan(stdev0), U = U*stdev0/stdev; Y = Kx+U; end
Z = U/stdev;

disp('Simulation setup finished')
