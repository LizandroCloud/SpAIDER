
function [xopt, fopt] = optfunc1Dinvecfibon(f,varargin)

% optfunc1Dinvecfibon -- ThreshLab's 1D optimization routine for 2D functions
%  Usage
%    [xopt fopt] = ...
%  optfunc1Dinvecfibon(f,'component',c,'starting',x0,'lo',xlo,'up',xup,...
%                      'eps',eps,par1,par2,...)
%  Inputs
%    f        : string of function name f(x;par1,par2,...)
%    x0       : starting point
%    xlo      : lower bound of search interval [xlo xup] (default 0)
%    xup      : upper bound of search interval [xlo xup] (default 1)
%    eps      : relative precision of output
%    component: index into x over which optimization takes place (default is 1)
%    'max' or 'min' : maximization or minimization (default is 'min')
%    par1, par1, ...: fixed parameters of f(x;par1,par2,...)
%  Outputs
%    xopt   : arg min/max f(x,par1,par2,...) for {x|x(c) in [xlo xup]}
%    fopt   : f(xopt)
%  Description
%    This is a simple (unconstrained) 1-d optimization procedure of the
%    function f(x;par1,par2,...), where par1, par2, ... are fixed parameter
%    values and x is an n-dimensional vector. The optimization takes place
%    over the c-th component of x, x(c), with 'component' c in
%    {1,...,length(x)}. All other components i are kept fixed, i.e., 
%    x(i) = x0(i)
%
%    The routine is based on Fibonacci optimization, it does not require
%    information on the derivative (gradient) but the object function (f)
%    should be (at least approximately) convex.
%  Note
%    Users of the optim(ization) toolbox may wish to replace the code by a
%    simple call of the functions fmin, fmincon or constr.
%  See also
%    help optfunc1Dfibon
%    help optfunc2Dinvecfibon

xlo = 0; xup = 1; eps = 1.e-4;
% x0 = Inf;
component = 1;
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
   case {'co'},
      component = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'ma'}
      minimize = 0;
      inputsinvarargin(k) = 0;
   case {'mi'}
      inputsinvarargin(k) = 0;
   end
end, end
if exist('x0') ~= 1, error('starting point must be provided'), end

if xup < xlo, xuplo = xup; xup = xlo; xlo = xuplo; end
% if x0 > xup, xup = 2*x0 - xlo; end
% if x0 < xlo, xlo = 2*x0 - xup; end


inputsinvarargin = inputsinvarargin(find(inputsinvarargin));
arginf = varargin(inputsinvarargin);

%initialisation Fibonacci numbers
F0 = 1;
F1 = 1;
F2 = 2;
Niter = 1;

% computation of number of iteration steps in optimisation
while (xup - xlo) / (xup + xlo) /F2 > eps
  F0 = F1;
  F1 = F2;
  F2 = F0 + F1;
  Niter = Niter + 1;
end

a = x0;
b = x0;
a(component) = xlo;
b(component) = xup;

v = a + (F1/F2) * (b - a);
fv = feval(f,v,arginf{:}); 
u = a + (F0/F2) * (b - a);
fu = feval(f,u,arginf{:});
if ~minimize, fu = -fu; fv = -fv; end
     
for k = 1:Niter - 1
   if (fu > fv)
      a = u;
      u = v;
      fu = fv;
      v = b-u+a;
      fv = feval(f,v,arginf{:});
      if ~minimize, fv = -fv; end
   else
      b = v;
      v = u;
      fv = fu;
      u = a+b-v;
      fu = feval(f,u,arginf{:});
      if ~minimize, fu = -fu; end
   end
end
     
xopt = (a+b)/2;
fopt = feval(f,xopt,arginf{:});
