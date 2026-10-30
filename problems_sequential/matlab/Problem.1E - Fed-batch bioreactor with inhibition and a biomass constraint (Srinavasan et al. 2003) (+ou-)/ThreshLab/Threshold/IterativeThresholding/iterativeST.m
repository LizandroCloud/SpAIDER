
function vhat = iterativeST(W,K,thr,maxit,vhat0);

% iterativeST -- Iterative soft thresholding
%  Usage
%    vhat = iterativeST(W,K,thr,maxit,vhat0);
%  Inputs
%    W      observations in sparsity model W = K*v+noise 
%    K      model of covariates (rectangular matrix)
%    thr    fixed threshold used in iterations
%    maxit  maximum number of iterations (default is 1000)
%    vhat0  initial value for vhat (default is zero)
%  Outputs
%    vhat   estimate of sparse vector v
%  Description
%    Iterative soft thresholding:
%    vhat_new = ST(vhat_old + K'*(W-K*vhat_old),thr)
%  Note
%  See also:
%    help iterativeHT
%    help iterativeGCVST
%    help iterativeMSEST
%    help LARS

wsize = size(W);
W = column(W);

if nargin < 2,
   K = eye(length(W));
end
[n m] = size(K);
if nargin < 3,
   thr = sqrt(2*log(m));
end

if nargin < 4,
   maxit = 1000;
end

if nargin < 5,
   vhat0 = zeros(m,1);
end

v0 = ones(m,1);
v1 = vhat0;

if norm(v1-v0) < eps, % i.e., if vhat0 happens to be ones(m,1)
   v0 = zeros(m,1);
end

lambda = thr;

i = 1;
while norm(v1-v0) > eps,
   v0 = v1;
   w1 = v0 + K'*(W-K*v0);
   v1 = ST(w1,lambda);
   i = i+1; 
   if i>maxit, 
      v0 = v1; 
      warning('iterativeST: maximum number of iterations reached')
   end
end
vhat = v1;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 

