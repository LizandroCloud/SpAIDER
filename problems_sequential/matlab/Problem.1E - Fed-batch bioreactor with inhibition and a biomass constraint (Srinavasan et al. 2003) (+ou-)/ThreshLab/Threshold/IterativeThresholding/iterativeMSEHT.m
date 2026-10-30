
function [vhat thr] = iterativeMSEHT(W,K,v,maxit,vhat0,KTK);

% iterativeMSEHT -- Iterative hard thresholding with step dependent threshold
%  Usage
%    [vhat thr] = iterativeMSEHT(W,K,v,maxit,vhat0,KTK);
%  Inputs
%    W      observations in sparsity model W = K*v+noise
%    K      model of covariates (rectangular matrix)
%    v      unobserved sparse vector v
%    maxit  maximum number of iterations
%    vhat0  initial value for vhat
%    KTK    binary var.(default 0); if one, then thresholding on standardized
%           components
%  Outputs
%    vhat   estimate of sparse vector v
%    thr    vector of threshold values used throughout iteration steps
%  Description
%    Iterative hard thresholding with step-dependent, minimum MSE threshold 
%    vhat_(i+1) = HT(vhat_(i) + K'*(W-K*vhat_old,thr(i))
%    or (when KTK == 1)
%    vhat_(i+1) = D*ST(inv(D)*((vhat_(i) + K'*(W-K*vhat_old,thr(i)))
%    with D = diag(sqrt(diag(K'*K)));
%  Note
%  See also:
%    help iterativeHT
%    help iterativeGCVHT


wsize = size(W);
W = column(W);

if nargin < 2,
   K = eye(length(W));
end
[n m] = size(K);

if nargin < 4,
   maxit = 1000;
end

if nargin < 5,
   vhat0 = zeros(m,1);
end

if nargin < 6,
   KTK = 0;
end

v0 = ones(m,1);
v1 = vhat0;

if norm(v1-v0) < eps, % i.e., if vhat0 happens to be ones(m,1)
   v0 = zeros(m,1);
end

D = ones(size(v0));
if KTK,
   D = sqrt(diag(K'*K));
end
Dinv = 1./D;
Dv = Dinv.*v;

i = 1;
thr = [];
while norm(v1-v0) > eps,
   v0 = v1;
   w1 = v0 + K'*(W-K*v0);
   w1 = Dinv.*w1;
   [v1 t] = minmsehard(Dv,w1);
   v1 = D.*v1;
   thr = [thr t];
   i = i+1; 
   if i>maxit, 
      v0 = v1; 
      warning('iterativeMSEHT: maximum number of iterations reached')
   end
end
vhat = v1;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
