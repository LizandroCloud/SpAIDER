
function [mse vhat] = MSEitthresh(v,w,K,t,maxit,vhat0,softORhard);

% MSEitthresh -- apply iterative soft/hard threshold 
%                + compute mean squared error (MSE)
%  Usage
%    [mse vhat] = MSEitthresh(v,K,w,t,softORhard,vhat0);
%  Inputs
%    v      noise-free sparse data
%    w      observed vector: w = K*v + noise
%    K      model of covariates (rectangular matrix)
%    t      threshold (must be scalar)
%    maxit  maximum number of iterations (default is 1000)
%    vhat0  initial value for vhat (default is zero)
%    softORhard boolean: 1 for soft; 0 for hard (default is 1)
%  Outputs
%    mse    the value of the mse of vhat
%    vhat   estimate of sparse vector v
%  Description
%    Applies soft- or hard threshold t to all elements of w and computes MSE
%  Note/See also
%    For minimizing MSEthresh, use minmseitthreshest

[n m] = size(K);
if nargin < 7, softORhard = 1; end
if nargin < 6, vhat0 = zeros(m,1); end
if nargin < 5, maxit = 1000; end

if isnan(maxit), maxit = 1000; end
if isnan(vhat0) | vhat0 == 0, vhat0 = zeros(m,1); end

if softORhard == 0,
   vhat = iterativeHT(w,K,t,maxit,vhat0);
else
   vhat = iterativeST(w,K,t,maxit,vhat0);
end
mse = norm(vhat-v)^2/m;
