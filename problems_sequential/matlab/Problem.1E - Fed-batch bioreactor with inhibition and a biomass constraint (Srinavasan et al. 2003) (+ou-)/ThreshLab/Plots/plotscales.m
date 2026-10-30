
function dummy = plotscales(w,L,wmin,wmax,varargin);

% plotscales - ThreshLab/Plots -- plot wavelet coefficients and borders scales
%   Usage
%     plotscales(w,L,wmin,wmax);
%     plotscales(w,L);
%     plotscales(w,L,wmin,wmax,option);
%     plotscales(w,L,option);
%   Inputs
%     w      sequence of wavelet coefficients
%     L      lowest scale
%     wmin   lower bound on y-axis of plot
%     wmax   upper bound on y-axis of plot
%     option plot options, as in matlab command plot
%   Outputs
%     void
%   Description
%     plots wavelet coefficients in w, with dashed lines indicating boundaries
%     between scales

wmindefined = 0;
if nargin == 1,
   options = {'b'};
   L = 0;
elseif nargin == 2, 
   if isa(L,'numeric'),
      options = {'b'};
   else
      options = {L};
      L = 0;
   end
elseif nargin == 3, 
   if isa(L,'numeric') & isa(wmin,'numeric'),
      wmax = wmin; wmin = L; L = 0; options = {'b'}; wmindefined = 1;
   elseif ~isa(L,'numeric'),
      options{1} = L; options{2} = wmin; L = 0;
   else
      options = {wmin};
   end
elseif nargin == 4,
   if ~isa(wmin,'numeric'),
      options{1} = wmin; options{2} = wmax;
   else
      options = {'b'};
      wmindefined = 1;
   end
else
   if isa(wmin,'numeric') & isa(wmax,'numeric'),
      for a = 1:nargin-4, options{a} = varargin{a}; end
      wmindefined = 1;
   else
      options{1} = wmin; options{2} = wmax;
      for a = 1:nargin-4, options{a+2} = varargin{a}; end
      wmindefined = 0;
   end
end
 
h = ishold;
N = length(w);
J = ceil(log2(N));
if ~wmindefined,
   if J == L,
      wmax = max(abs(w));
      wmin = -wmax;
   else
      wmax = max(abs(w(2^L+1:N)));
      wmin = -wmax;
   end
end
vert = [wmin,wmax];
plotasspikes(w,options{:})
hold on
scaleborder = [N N];
for j = J-1:-1:L
   scaleborder = scaleborder/2;
   plot(scaleborder,vert,'--')
end
axis([0,N,wmin,wmax])
if ~h,
   hold off  
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
