
function [xopt fopt] = optfunc1Dfibon(f,varargin)

% optfunc1Dfibon -- ThreshLab's optimization routine for 1D functions
%  Usage
%    [xopt fopt] = ...
%  optfunc1Dfibon(f,'lo',xlo,'up',xup,'eps',eps,'argpos',argpos,par1,par2,...)
%  Inputs
%    f      : string of function name
%    xlo    : lower bound of search interval [xlo xup] (default 0)
%    xup    : upper bound of search interval [xlo xup] (default 1)
%    eps    : relative precision of output
%    argpos : position of variable in parameter list (default is 1)
%    'max' or 'min' : maximization or minimization (default is 'min')
%    par1, par1, ...: fixed parameters of f(x;par1,par2,...)
%  Outputs
%    xopt   : arg min/max f(x) for x in [xlo xup]
%    fopt   : f(xopt)
%  Description
%    This is a simple (unconstrained) 1-d optimization procedure of the
%    function f(x;par1,par2,...), where par1, par2, ... are fixed parameter
%    values and x is a 1-d variable over which the optimization takes place.
%
%    if argpos>1, then the function can be optimized w.r.t. other than its
%    first argument. E.g., if argpos = 3, then the routine optimizes
%    f(par1,par2,x,par3,...). Argpos must be smaller than input arguments of f
% 
%    The routine is based on Fibonacci optimization, it does not require
%    information on the derivative (gradient) but the object function (f)
%    should be (at least approximately) convex.
%  Note
%    Users of the optim(ization) toolbox may wish to replace the code by a
%    simple call of the functions fmin, fmincon or constr.
%  See also
%    help minfunc1Dfibon
%    help maxfunc1Dfibon

xlo = 0; xup = 1; eps = 1.e-4;
% x0 = Inf;
argpos = 1;
nvarargin = length(varargin);
inputsinvarargin = (1:nvarargin);
minimize = 1;
for k = 1:nvarargin, if ischar(varargin{k}), 
   varname = varargin{k};
   switch varname(1:2),
   case {'lo','xl'},
      xlo = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'hi','xh','up','xu'},
      xup = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'ep'},
      eps = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'x0','st','in'},
      x0 = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'ar'},
      argpos = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'ma'}
      minimize = 0;
      inputsinvarargin(k) = 0;
   case {'mi'}
      inputsinvarargin(k) = 0;
   end
end, end
% if isinf(x0), x0 = (xlo+xup)/2; end

if xup < xlo, xuplo = xup; xup = xlo; xlo = xuplo; end
% if x0 > xup, xup = 2*x0 - xlo; end
% if x0 < xlo, xlo = 2*x0 - xup; end

% x0 is not used in this version (can be omitted)

inputsinvarargin = inputsinvarargin(find(inputsinvarargin));
parinf = varargin(inputsinvarargin);
narginf = length(parinf)+1;
arginf = cell(1,narginf);
for i = 1:argpos-1, arginf{i} = parinf{i}; end
for i = argpos+1:narginf, arginf{i} = parinf{i-1}; end

%initialisation Fibonacci numbers
F0 = 1;
F1 = 1;
F2 = 2;
Niter = 1;

% computation of number of iteration steps in optimisation
while (xup - xlo) / F2 > eps
  F0 = F1;
  F1 = F2;
  F2 = F0 + F1;
  Niter = Niter + 1;
end


% optimization

a = xlo;
b = xup;

v = a + (F1/F2) * (b - a); arginf{argpos} = v;
fv = feval(f,arginf{:});
u = a + (F0/F2) * (b - a); arginf{argpos} = u;
fu = feval(f,arginf{:});
if ~minimize, fu = -fu; fv = -fv; end

for k = 1:Niter - 1
   if (fu > fv)
      a = u;
      u = v;
      fu = fv;
      v = b-u+a; arginf{argpos} = v;
      fv = feval(f,arginf{:});
      if ~minimize, fv = -fv; end
   else
      b = v;
      v = u;
      fv = fu;
      u = a+b-v; arginf{argpos} = u;
      fu = feval(f,arginf{:});
      if ~minimize, fu = -fu; end
   end
end

xopt = (a+b)/2; arginf{argpos} = xopt;
fopt = feval(f,arginf{:});

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 

