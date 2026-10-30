function s = shift(v,l)

% shift - ThreshLab/General/ -- circular shift of matrix or vector
%   Usage
%     s = shift(v,l);
%   Inputs
%     v   : vector or matrix
%     l   : integer number (default 1)
%   Outputs
%     s   : shifted version of v
%   Description
%     shifts columns of matrix v or elements of vector v by l positions to the
%     left
%     l may be positive or negative

if nargin == 1,
   l = 1;
end
vv = v;
if size(v,2) == 1, v = v'; end

l = mod(l,size(v,2));

if isempty(v),
   s = v;
else
   s = [v(:,l+1:size(v,2)) v(:,1:l)];
   s = reshape(s,size(vv));
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
