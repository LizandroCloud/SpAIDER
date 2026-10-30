
function [xopt, fopt] = optfunc2Dinvecfibon(f,varargin)

% optfunc2Dinvecfibon -- ThreshLab's optimization routine for 2D or nD functions
%  Usage
%    [xopt fopt] = ...
%  optfunc2Dinvecfibon(f,'starting',x0,'lo',xlo,'up',xup,...
%                 'eps',eps,par1,par2,...)
%  Inputs
%    f        : string of function name f(x;par1,par2,...)
%    x0       : starting point
%    xlo      : lower bound of search interval [xlo xup] (default (0,0,0,...))
%    xup      : upper bound of search interval [xlo xup] (default (1,1,1,...))
%    eps      : relative precision of output
%    'max' or 'min' : maximization or minimization (default is 'min')
%    par1, par1, ...: fixed parameters of f(x;par1,par2,...)
%  Outputs
%    xopt   : arg min/max f(x,par1,par2,...) for {x in [xlo xup]}
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
%    help optfunc1Dinvecfibon


eps = 1.e-4;
nvarargin = length(varargin);
inputsinvarargin = (1:nvarargin);
for k = 1:nvarargin, if ischar(varargin{k}),
   varname = varargin{k};
   switch varname(1:2),
   case {'lo','xl'},
      xlo = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'hi','xh','up','xu'},
      xup = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'x0','st','in'},
      x0 = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   end
end, end
if exist('x0') ~= 1, error('starting point must be provided'), end

inputsinvarargin = inputsinvarargin(find(inputsinvarargin));
varargin1D = varargin(inputsinvarargin);


dim = length(x0);
if exist('xlo') ~= 1, xlo = zeros(1,dim); end
if exist('xup') ~= 1, xup = ones(1,dim); end

relerr = 1;
comp = 1;
while relerr > eps,
   [xopt fopt] = ...
   optfunc1Dinvecfibon(f,'lo',xlo(comp),'hi',xup(comp),'component',comp,...
                       'starting',x0,varargin1D{:});
   relerr = norm(xopt - x0)/(norm(xopt) + norm(x0));
   comp = mod(comp,dim)+1;
   x0 = xopt;
end

