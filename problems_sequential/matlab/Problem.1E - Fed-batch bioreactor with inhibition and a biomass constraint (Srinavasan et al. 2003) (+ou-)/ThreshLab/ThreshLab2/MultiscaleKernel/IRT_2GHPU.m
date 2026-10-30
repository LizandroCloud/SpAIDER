
function y = IRT_2GHPU(x,w,M,varargin);

% RWT_2GHtPU -- ThreshLab2/Irregular/MultiscaleKernel -- ThreshLab's routine
%               for signal reconstruction from:
%                Burt-Adelson's pyramid on irregular point sets
%                Smoothed lifted smoothing 
%                Redundant Wavelet Transform (RWT), 2nd Generation (2G) with:
%                prefilter H, prediction P, update U
%  Usage
%    y = IRT_2GHPU(x,w,M,varargin);
%  Inputs
%    x      grid-locations of observations (time, space)
%    w      wavelet/pyramid coefficients
%    M      number of resulution levels
%    varargin several options (see below)
%  Outputs
%    y      vector of size(x) with reconstructed observations/signal
%    Ne     (optional) number of scaling coefficients (dependent on M)
%  Description
%    Thanks to the overcompleteness of the decomposition, this reconstruction
%    does not depend on the prefilter step in the forward transform.
%    The reconstruction uses (= undoes) one update step and prediction step.
%    These steps may be preceded by a prefilter step (which may be identical to
%    the prefilter in the decomposition)
%  Examples
%    Set up
%    Analyses
%         
%  See also
%    help IRT_2GHPU1step
%    help RWT_2GHtPU

x = column(x);
N = length(x);
ytranspose = 0;
if size(w,1) == 1, w = w'; ytranspose = 1; end
ncols = size(w,2);

n = 0; % default value; may be overwritten by subsequent statement
moments = [];
Ne = NaN;
nvarargin = length(varargin);
for k = 1:nvarargin-1, vark = varargin{k}; if ischar(vark), switch(vark)
    case {'n','nU','primom','primalmoments','nprimal','nprimalmoments'}
       n = varargin{k+1};
    case {'moments','Moments'}
       moments = varargin{k+1};
    case {'Ne'}
       Ne = varargin{k+1};
end, end, end

if isempty(moments) & n > 0,
   [void, Ne, moments] = ...
   RWT_2GHtPU(x,sparse(zeros(size(x))),M,varargin{:},'U',0,'Ht',1);
   % options in varargin w.r.t. U and Ht are overwritten by adding U=0 and Ht=1
   % As only P determines the moments, we save computations with U=0 and Ht=1.
end

% Now, find grids at all levels
even = (1:N);
gridj = even;
Nje = zeros(1,M+1);
Nje(M+1) = N;
for j=M:-1:1,
   even = even(splitevenodd(even));
   gridj = [even gridj];
   Nje(j) = length(even);
end

onestepvarargin = varargin;
% every step in the reconstruction overwrites the details and throws away the
% evens (see flow diagrams)

Ne = Nje(1);
coarse = gridj(1:Ne); 
xcoarse = x(coarse);
if n>0,
   Mcoarse = moments(1:Ne,1:n);
else
   Mcoarse = [];
end
sj = w(1:Ne,1:ncols);
k0 = Ne;

for j = 1:M,
   No = Nje(j+1);
   fine = gridj(k0+1:k0+No);
   even = splitevenodd(fine);
   xfine = x(fine);
   if n>0,
      Mfine = moments(k0+1:k0+No,1:n);
   else
      Mfine = [];
   end
   wj = w(k0+1:k0+No,1:ncols);
   bandwidthfactor = 2^(M-j);
   [sj varargout] = IRT_2GHPU1step(xfine,even,wj,sj,Mfine,Mcoarse,...
                       onestepvarargin{:},...
                       'bandwidthfactor',bandwidthfactor);
   % prepare for next step
   Ne = No; k0 = k0+No;
   coarse = fine; xcoarse = xfine; Mcoarse = Mfine;
   onestepvarargin = varargout;
end

y = sj;
if ytranspose, y = y'; end
