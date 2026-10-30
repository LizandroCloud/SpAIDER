
function [sfine listofoptions] = ...
          IRT_2GHPU1step(xfine,even,wcoarse,scoarse,Mfine,Mcoarse,varargin);

% IRT_2GHPU1step -- ThreshLab2/Irregular/MultiscaleKernel -- One step in
%                ThreshLab's routine for signal reconstruction from:
%                Burt-Adelson's pyramid on irregular point sets
%                Smoothed lifted smoothing 
%                Redundant Wavelet Transform (RWT), 2nd Generation (2G) with:
%                prefilter H, prediction P, update U
%  Usage
%     [sfine listofoptions] =
%     IRT_2GHPU1step(xfine,xcoarse,wj,sj,Mfine,Mcoarse,varargin);
%  Inputs
%    xfine    fine scale grid-locations of observations (time, space)
%    even     index vector into xfine for coarse scale grid-locations
%    wcoarse  detail/wavelet coefficients or offsets at coarse scale
%    scoarse  scaling coefficients at coarse scale
%    Mfine    fine scale scaling function moments
%    Mcoarse  coarse scale scaling function moments
%  Outputs
%    sfine    fine scale scaling coefficients
%  Description
%    Applies one step inverse lifting filter bank with an optional prefilter,
%    at most one prediction step and at most one update step. 
%  See also
%    help IRT_2GHPU

if size(xfine,2) > size(xfine,1),
   xfine = xfine';
end
dotranspose = 0;
if size(wcoarse,1) == size(scoarse,1),
   wcoarse = wcoarse';
   scoarse = scoarse';
   dotranspose = 1;
end

x = xfine; xeven = x(even);
N = length(x);
n = size(Mfine,2);
% following values are the default values; they may be overwritten while
% parsing through varargin
H = 0;
P = 1;
U = 0;
nt = 1;
bandwidths = (max(x) - min(x)) * 3/N;
bandwidthfactor = 1;
kerneltype = 'cos';
prefilteroptions = {};

% Parsing varargin
nvarargin = length(varargin);
listofoptions = {};
for k = 1:nvarargin-1, vark = varargin{k}; if ischar(vark), switch(vark)
   case {'n','nU','primom','primalmoments','nprimal','nprimalmoments'}
      if n ~= varargin{k+1},
   warning('Inconsistent # prim. moments ignored; size(Mfine,2) used instead')
      end
   case {'nt','ntilde','ntP','dualmom','dualmoments','ndual','ndualmoments'}
      nt = varargin{k+1};
      listofoptions = {listofoptions{:},varargin{k:k+1}};
   case {'H','Hfilter','prefilter'}
      H = varargin{k+1};
      listofoptions = {listofoptions{:},varargin{k:k+1}};
   case {'P','Pmethod','prediction','Prediction','Predict','predict'}
      P = varargin{k+1};
      listofoptions = {listofoptions{:},varargin{k:k+1}};
   case {'U','update'}
      U = varargin{k+1};
      listofoptions = {listofoptions{:},varargin{k:k+1}};
   case {'bandwidth','bandwidths','h','hj','hJ'}
      bandwidths = varargin{k+1};
      % bandwidth in listofoptions will be different (see below)
   case {'bandwidthfactor'}
      bandwidthfactor = varargin{k+1};
      % bandwidthfacor will not be passed tothe next scale (calling routine has
      % to provide a new value)
   case {'kerneltype','type'}
      kerneltype = varargin{k+1};
      listofoptions = {listofoptions{:},varargin{k:k+1}};
   case {'prefilteroptions','prefiltoptions','Htopt','prefiltoption'}
      Htoption = varargin{k+1};
      if iscell(Htoption),
         prefilteroptions = {prefilteroptions{:},Htoption{:}};
      else
         prefilteroptions = {prefilteroptions{:},Htoption};
      end
      listofoptions = {listofoptions{:},varargin{k:k+1}};
end, end, end

% Design of the transform
  % PREDICTION MATRIX %
  Ne = length(xeven);
  Ie = sparse(eye(Ne));
  if ischar(P), switch(P(1:6))
  case {'kernel'}
     nh = length(bandwidths);
     h =  bandwidths(nh)*bandwidthfactor;
     P = sparse(kernelestimation01(x,xeven,Ie,h,kerneltype)');
  case {'inter','polyno','intpol'}
     P = sparse(polynomialprediction(xeven,Ie,x,nt));
  case {'locpol','localp'}
     nh = length(bandwidths);
     h =  bandwidths(nh)*bandwidthfactor;
     P = sparse(locpolestimation00(x,xeven,Ie,h,kerneltype,nt)');
  otherwise
     P = 1;
  end, end
  if isscalar(P), P = sparse(diag(P*ones(1,N))); P = P(1:N,even); end
  % PREFILTER %
  if ischar(H),
     if H(2:5) == 'iort',
        H = findbiorthfilter(P,Mfine);
     else
        H = constructprefilter(H,xfine,prefilteroptions);
     end
  end
  % UPDATE %
  if n>0,
     Mf = Mfine(1:N,1:n); Mc = Mcoarse(1:Ne,1:n);
     Io = sparse(eye(N));
     U = updateprimalmoments(x,Io,xeven,Mf,Mc);
  else
     U = 0;
  end
  if isscalar(U), U = sparse(diag(U*ones(1,N))); U = U(even,1:N); end
% Execution of transforms
  % PREFILT %
  Hwcoarse = H*wcoarse;
  % SUBSAMPLE AND SUBTRACT%
  scoarse = scoarse-Hwcoarse(even,:);
  % UPDATE %
  scoarse = scoarse - U*wcoarse;
  % PREDICT %
  wcoarse = wcoarse + P*scoarse;

% Filling in
% Prepare next step
  nh = length(bandwidths);
  if nh > 1, 
     bandwidths = bandwidths(1:nh-1);
  end
  listofoptions = {listofoptions{:},'bandwidths',bandwidths};
  sfine = wcoarse;
