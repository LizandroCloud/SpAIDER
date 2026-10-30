
function Q = gramschmidt(A)

% ThreshLab/General/gramschmidt.m -- ThreshLab's Gram Schmidt orthogonalization
%  Usage
%    Q = gramschmidt(A)
%  Input
%    A   matrix
%  Output
%    Q   matrix with same size as A
%   Description
%    if size(A) = [n m], then the first min(m,n) columns of Q are orthogonal 

Q = ones(size(A));
m = size(A,2);
Q(:,1) = A(:,1)/norm(A(:,1));
for k = 2:m
  c = Q(:,1:k-1)'*A(:,k);
  V = A(:,k) - Q(:,1:k-1)*c;
  if norm(V) < eps,
     Q(:,k) = 0;
  else
     Q(:,k) = V/norm(V);
  end
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
