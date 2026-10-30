
function Y = metropolisIsing(X,W,sig,tau,thr)

% metropolisIsing -- one step of metropolis sampling for Ising Random Field
%   Usage
%    Y = metropolisIsing(X,W,sig,tau,thr);
%   Inputs
%    X    state labels in Markov chain Monte Carlo (MCMC) iteration step
%    W    observations
%    sig  standard deviation of observational noise
%    tau  parameter of Ising model (see below)
%    thr  threshold
%   Outputs
%    Y    state labels in next MCMC iteration step
%   Description
%    Y is a random perturbation of X, taking into account the values of W
%    We observe W = V + noise, where V has zero-inflated maixture distribution
%    X is the hidden/latent binary state vector of V (zero or non zero)
%    model for X is Ising model, with "temperature" tau
%    threshold thr is used to estimate X from W: X = HT(W,thr)
%   Note
%    This is a vector implementation, which means that update of one component
%    is not yet active in updating other components in this iteration step
%    (slower convergence rate per iteration, but faster iterations)
%   See also
%    help isingpotupdate
%    help condprob2

[N M] = size(X);
X1 = 1 - X;
X0 = (abs(W)>thr);
p = sum(sum(X0))/(N*M);
eps = 1.e-9;
p = min(p,1-eps);
p = max(p,eps);
condratio = condprob2(1-X,W,sig,thr,p)./condprob2(X,W,sig,thr,p);
priordif = isingpotupdate(2*X-1,1-2*X,tau) - isingpotupdate(2*X-1,2*X-1,tau);
XX = X(2:N-1,2:M-1);
r = exp(-priordif).*condratio(2:N-1,2:M-1);
Y = X;
Y(2:N-1,2:M-1) = (rand(size(r)) <= r).*(1 - 2*X(2:N-1,2:M-1)) + X(2:N-1,2:M-1);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
