
function [betahat, muhat, x, Cp, lambda] = LARS(Y,K,stdev,n1)

% LARS - Least Angle Regression for LASSO (beta-version; subject to change)
%  Usage
%    [betahat muhat x Cp] = LARS(Y,K,stdev,n1)
%  Inputs
%    Y      observations in sparsity model Y = K*beta+noise
%    K      model of covariates (rectangular matrix)
%    stdev  (estimated) standard deviation (default is 1; current version does
%           not estimate stdev)
%    n1     desired number of nonzeros (default is n=length(Y))
%  Outputs
%    betahat estimator beta
%    muhat   estimator mu = K*beta
%    x       ordered set of selected variables (indices in {1,..,length(beta)})
%    Cp      Mallows' Cp values of nested selections
%    lambda  values of penalty parameter corresponding to x and Cp
%  Description
%  Note
%     see p.74 in Elements of statistical learning, second edition;
%     downloadable from http://www-stat.stanford.edu/~tibs/ElemStatLearn/
%     see also paper LARS, Ann.Stat.2004, pages 407-499
%  See also:
%    help iterativeST

if nargin < 4, n1 = n; end
if nargin < 3, stdev = 1; end

Y = column(Y);

ybar = mean(Y);
Y = Y-ybar;
muhat = zeros(size(Y));
n = length(Y);
m = size(K,2);
colnormK = sqrt(sum(K.^2,1));
X = K*diag(1./colnormK);

i = 0;
dof = i;
% Cp = [Cp(-1) Cp(0) Cp(1)]
Cp = [Inf norm(Y-muhat)^2/n+2*dof/n*stdev^2-stdev^2];
lambda = [Inf Inf];
A = [];
Aprime = (1:m);
% while Cp(i+1) > Cp(i+2),
while i < n1,
    muhat0 = muhat;
    chat = X'*(Y-muhat0);
    [Chat j] = max(abs(chat));
    if i>0,
       s = sign(chat(A));
       XA = X(1:n,A)*diag(s);
       GA = XA'*XA;
       AA = 1/sqrt(ones(1,i)*(GA\ones(i,1)));
       wA = AA*(GA\ones(i,1));
       uA = XA*wA;
    else
       uA = X(1:n,j); uA = uA/norm(uA);
       AA = Inf;
    end
    a = X'*uA;
    gamma = [(Chat-chat(Aprime))./(AA-a(Aprime)) ...
             (Chat+chat(Aprime))./(AA+a(Aprime))];
    gamma(gamma<0) = Inf;
    gamma = min(gamma,[],2);
    [gammahat j] = min(gamma);
    muhat = muhat0 + gammahat*uA;
    newAprime = [Aprime(1:j-1) Aprime(j+1:m-i)]; 
    j = Aprime(j); A = [A j]; Aprime = newAprime;
    i = i+1;
    dof = i;
    lambda = [lambda Chat];
    Chat = Chat -gammahat*AA;
    Cp = [Cp norm(Y-muhat0)^2/n+2*dof/n*stdev^2-stdev^2];
end
muhat = muhat0+ybar;
A = A(1:i);
XA = X(1:n,A);

% HT
% betahatA = (XA'*XA)\(XA'*Y);
% ST
betahatA = (XA'*XA)\(XA'*muhat);
betahat = zeros(1,m);
betahat(A) = betahatA;

x = A;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
