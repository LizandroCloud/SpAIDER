
function [xmax fmax] = maxfunc1Dfibon(f,varargin)

% maxfunc1Dfibon -- ThreshLab's maximization routine for 1D functions
%  Usage
%    [xmax fmax] = ...
%  maxfunc1Dfibon(f,'lo',xlo,'up',xup,'eps',eps,'argpos',argpos,par1,par2,...)
%  Inputs
%    f      : string of function name
%    xlo    : lower bound of search interval [xlo xup] (default 0)
%    xup    : upper bound of search interval [xlo xup] (default 1)
%    eps    : relative precision of output
%    argpos : position of variable in parameter list (default is 1)
%    par1, par1, ...: fixed parameters of f(x;par1,par2,...)
%  Outputs
%    xmax   : arg max f(x) for x in [xlo xup]
%    fmax   : f(xmax)
%  Description
%    This is a simple (unconstrained) 1-d maximization procedure of the
%    function f(x;par1,par2,...), where par1, par2, ... are fixed parameter
%    values and x is a 1-d variable over which the optimization takes place.
%
%    if argpos>1, then the function can be maximized w.r.t. other than its
%    first argument. E.g., if argpos = 3, then the routine maximizes
%    f(par1,par2,x,par3,...). Argpos must be smaller than input arguments of f
% 
%    The routine is based on Fibonacci maximization, it does not require
%    information on the derivative (gradient) but the object function (f)
%    should be (at least approximately) convex.
%  Note
%    Users of the optim(ization) toolbox may wish to replace the code by a
%    simple call of the functions fmax, fmaxcon or constr.
%  See also

varargin2 = {varargin{:} 'max'};
[xmax fmax] = optfunc1Dfibon(f,varargin2{:});
