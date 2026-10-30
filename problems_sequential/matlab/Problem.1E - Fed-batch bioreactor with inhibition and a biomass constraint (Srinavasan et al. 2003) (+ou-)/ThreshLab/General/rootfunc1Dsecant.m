
function x = rootfunc1Dsecant(f,varargin)

% rootfunc1Dsecant -- ThreshLab's secant method for finding roots of 1D func.
%  Usage
%    x = ...
%  rootfunc1Dsecant(f,'x0',x0,'x1',x1,'eps',eps,'maxiter',maxit,'argpos',argpos,'fval',fo,par1,par2,...)
%  Inputs
%    f      : string of function name
%    x0     : initial value (default 0)
%    x1     : second initial value (default 1)
%    eps    : precision of output
%    maxit  : maximum number of iterations (default 1000)
%    argpos : position of variable in parameter list (default is 1)
%    par1, par1, ...: fixed parameters of f(x;par1,par2,...)
%  Outputs
%    x      : zero/root of f
%  Description
%  Note
%    Users of the funfun toolbox may wish to replace this routine by a call of
%    fzero
%  See also
%    help optfunc1Dfibon

x0 = 0; x1 = 1;
fo = 0;
eps = 1.e-4;
maxit = 1000;
argpos = 1;

nvarargin = length(varargin);
inputsinvarargin = (1:nvarargin);
minimize = 1;
for k = 1:nvarargin, if ischar(varargin{k}), 
   varname = varargin{k};
   switch varname(1:2),
   case {'lo','xl','x0'},
      x0 = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'hi','xh','up','xu','x1'},
      x1 = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'ep'},
      eps = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'ar'},
      argpos = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'fv','fo'},
      fo = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   case {'ma'},
      maxit = varargin{k+1};
      inputsinvarargin(k:k+1) = 0;
   end
end, end

inputsinvarargin = inputsinvarargin(find(inputsinvarargin));
parinf = varargin(inputsinvarargin);
narginf = length(parinf)+1;
arginf = cell(1,narginf);
for i = 1:argpos-1, arginf{i} = parinf{i}; end
for i = argpos+1:narginf, arginf{i} = parinf{i-1}; end


% iteration
iter = 0;

arginf{argpos} = x0;
f0 = feval(f,arginf{:});
arginf{argpos} = x1;
f1 = feval(f,arginf{:});
xx = x0;
ff = f0;
if abs(f1-fo) < abs(ff-fo), xx = x1; ff = f1; end

while abs(ff-fo) > eps,
   iter = iter+1;
   if iter > maxit,
      ff = fo;
      xx = NaN;
   else
      x2 = x1 - (f1-fo)*(x1-x0)/(f1-f0);
      arginf{argpos} = x2;
      f2 = feval(f,arginf{:});
      if abs(f2-fo) < abs(ff-fo), xx = x2; ff = f2; end
      f0 = f1; x0 = x1; f1 = f2; x1 = x2;
   end
end

x = xx;
