
function [w, Ne, moments] = RWT_2GHtPU(x,y,M,varargin);

% RWT_2GHtPU -- ThreshLab2/Irregular/MultiscaleKernel -- ThreshLab's routine 
%               for:
%                Burt-Adelson's pyramid on irregular point sets
%                Smoothed lifted smoothing 
%                Redundant Wavelet Transform (RWT), 2nd Generation (2G) with:
%                prefilter Ht, prediction P, update U
%  Usage
%    [w, Ne, moments] = RWT_2GHtPU(x,y,M,varargin);
%  Inputs
%    x        grid-locations of observations (time, space)
%    y        observed values
%    M        number of resulution levels
%    varargin several options (see below)
%  Outputs
%    w        vector with wavelet coefficients
%    Ne       (optional) number of scaling coefficients (dependent on M)
%    moments  (optional) moments of intermediate scaling functions (before
%                       update in each step)
%  Description
%    Applies repeated smoothed lifting filter bank with prefiltering,
%    one prediction step and one update step. 
%  Examples
%    Set up
%       N = 1000;
%       x = makegrid(N,1);
%       f = makesignal(x,'skyline');
%       y = makenoisy(f,10);
%       n = 2; M = 7;
%       predkernel = 'cos';
%       finescalebandwidth = (max(x) - min(x)) * 3/N;
%       hJ = finescalebandwidth;
%    Analyses
%         
%  See also
%    help RWT_2GHtPU1step
%    help IRT_2GHPU

x = column(x);
N = length(x);
wtranspose = 0;
if size(y,1) ~= N, y = y'; wtranspose = 1; end
if size(y,1) ~= N, error('dimension mismatch x and y'), end
ncols = size(y,2);

n = 0; % default value; may be overwritten by subsequent statement
nvarargin = length(varargin);
for k = 1:nvarargin-1, vark = varargin{k}; if ischar(vark), switch(vark)
    case {'n','nU','primom','primalmoments','nprimal','nprimalmoments'}
       n = varargin{k+1};
end, end, end

if n>0,
   Mcoarse = initmom(x,(0:n-1));
else
   Mcoarse = [];
end
moments = Mcoarse;

xcoarse = x;
swcoarse = y;
w = swcoarse;
Ne = length(xcoarse);
Nw = size(w,1);
onestepvarargin = varargin;
for j = M:-1:1,
   xfine = xcoarse;
   sfine = swcoarse(1:Ne,1:ncols);
   Mfine = Mcoarse;
   previousdetails = w(Ne+1:Nw,1:ncols);
   [xcoarse swcoarse Mcoarse listofoptions] = ...
            RWT_2GHtPU1step(xfine,sfine,Mfine,onestepvarargin{:});
   w = [swcoarse; previousdetails];
   onestepvarargin = listofoptions;
   Nw = size(w,1);
   Ne = length(xcoarse);
   moments = [Mcoarse;moments];
end

if wtranspose, w = w'; end
