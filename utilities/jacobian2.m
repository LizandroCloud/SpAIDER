function J = jacobian(F,x)
%JACOBIAN Jacobian matrix.
%   JACOBIAN(f,v) computes the Jacobian of the scalar or vector f 
%   with respect to the vector v. The (i,j)-th entry of the result
%   is df(i)/dv(j). Note that when f is scalar, the Jacobian of f
%   is the gradient of f. Also, note that scalar v is allowed, 
%   although this is just DIFF(f,v).
%
%   Example:
%       jacobian([x*y*z; y; x+z],[x y z])
%       jacobian(u*exp(v),[u;v])
%
%   See also SYM/DIFF.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.10.4.3 $  $Date: 2006/06/20 20:52:19 $

if nargin == 1
   x = ['[' findsym(F) ']'];
end
J = sym('jacobian',rowv(F),rowv(x));


function v = rowv(x)
% Convert x to character row vector, eliminating any 'array([' and ']).
v = sym(x);
if all(size(v) == 1)
   v = ['[' char(v) ']'];
else
  % v = char(v,1);
   v = char(v);
   v(1:7) = [];
   v(end) = [];
end
