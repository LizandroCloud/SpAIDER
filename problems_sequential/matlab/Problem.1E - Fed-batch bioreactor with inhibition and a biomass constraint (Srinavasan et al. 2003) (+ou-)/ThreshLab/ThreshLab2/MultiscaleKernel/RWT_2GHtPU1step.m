
function [xcoarse wcoarse Mcoarse listofoptions] = ...
            RWT_2GHtPU1step(xfine,sfine,Mfine,varargin);

% RWT_2GHtPU1step -- ThreshLab2/Irregular/MultiscaleKernel -- One step in
%                ThreshLab's routine 
%                Burt-Adelson's pyramid on irregular point sets
%                Smoothed lifted smoothing
%                Redundant Wavelet Transform (RWT), 2nd Generation (2G) with:
%                prefilter Ht, prediction P, update U
%  Usage
%    [xcoarse wcoarse Mcoarse listofoptions] =
%    RWT_2GHtPU1step(xfine,sfine,Mfine,varargin);
%  Inputs
%    xfine  fine scale grid-locations of observations (time, space)
%    sfine  fine scale scaling coefficients
%    Mfine  fine scale moments of scaling basis functions
%    ...    several optional parameters, for use in constructprefilter
%    'n', n   n being the number of vanishing moments in the update step
%                  (primal moments); default value is size(Mfine,2)
%
%  Outputs
%    xcoarse coarse scale grid-locations
%    wcoarse coarse scale scaling and wavelet coefficients 
%    Mcoarse coarse scale moments of scaling basis functions
%  Description
%    Applies one step lifting filter bank with an optional prefilter, at most 
%    one prediction step and at most one update step. 
%  Note
%    the first length(xcoarse) elements in wcoarse are coarse scale scaling
%    coefficients. Remaining elements in wcoarse are coarse scale wavelet
%    coefficients
%  See also
%    help RWT_2GHtPU

transposex = 0;
if size(xfine,2) > size(xfine,1),
   xfine = xfine';
end
transposew = 0;
if size(sfine,1) ~= size(xfine,1),
   sfine = sfine';
end

x = xfine;
y = sfine;
N = length(x);
n = size(Mfine,2);
% following values are the default values; they may be overwritten while
% parsing through varargin
Ht = 1;
P = 1;
U = 0;
nt = 1;
bandwidths = (max(x) - min(x)) * 3/N;
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
   case {'Ht','Htfilter','prefilter'}
      Ht = varargin{k+1};
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
  % SPLIT %
  [even, odd] = splitevenodd(x);
  yodd  = y(odd,:);
  yeven = y(even,:);
  xodd = x(odd);
  xeven = x(even);
  % PREDICTION MATRIX %
  Ne = length(even);
  Ie = sparse(eye(Ne));
  if ischar(P), switch(P(1:6))
  case {'kernel'}
     h =  bandwidths(1);
     P = sparse(kernelestimation01(x,xeven,Ie,h,kerneltype)');
  case {'inter','polyno','intpol'}
     P = sparse(polynomialprediction(xeven,Ie,x,nt));
  case {'locpol','localp'}
     h =  bandwidths(1);
     P = sparse(locpolestimation00(x,xeven,Ie,h,kerneltype,nt)');
  otherwise
     P = 1;
  end, end
  if isscalar(P), P = sparse(diag(P*ones(1,N))); P = P(1:N,even); end
  % MOMENTS %
  % Mcoarse = Mfine(even,:) + P'*Mfine;
  if ~isempty(Mfine),
     Mcoarse = P'*Mfine;
  else
     Mcoarse = Mfine;
  end
  % PREFILTER %
  if ischar(Ht), Ht2 = [Ht '        '];
     if Ht2(2:5) == 'iort',
        Ht = findbiorthfilter(P,Mfine);
     else
        Ht = constructprefilter(Ht,xfine,prefilteroptions);
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
% Execution of transform
  % PREFILT %
  yHt = Ht*y;
  % SUBSAMPLE %
  yeven = yHt(even,:);
  % PREDICT %
  offsets = y - P*yeven;
  % UPDATE %
  yeven = yeven + U*offsets;

% Filling in
% Prepare next step
  nh = length(bandwidths);
  if nh == 1, bandwidths = bandwidths*2;
  else bandwidths = bandwidths(2:nh);
  end
  listofoptions = {listofoptions{:},'bandwidths',bandwidths};

  No = length(xodd);
  Ne = length(xeven);
  wcoarse = zeros(2*Ne+No,size(sfine,2));
  wcoarse(Ne+1:2*Ne+No,:) = offsets;
  wcoarse(1:Ne,:)         = yeven;

  xcoarse = xeven;

if transposex, xcoarse = xcoarse'; end
if transposew, wcoarse = wcoarse'; end
