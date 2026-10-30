
function H = constructprefilter(H,x,varargin);

% constructprefilter -- construct prefilter operation in smoothed lifting
%  Usage
%    H = constructprefilter(H,x);
%    H = constructprefilter(H,n);
%    H = constructprefilter(H,x,optional parameters);
%    H = constructprefilter(H,n,optional parameters);
%  Inputs
%    H   character string with name of filter operation
%    x   grid/locations of (mostly nonequidistant) observations
%    n   number of (equidistant) observations
%    optional parameters, depending on choice of H
%  Outputs
%    H   sparse filter matrix
%  Description
%  Examples
%    H = constructprefilter('I',x);
%        yields identical matrix
%    H = constructprefilter('Haar',x);
%        yields Haar smoothing (1/2,1/2)
%    H = constructprefilter('Unbalanced Haar',x,Mfine,Mcoarse);
%    H = constructprefilter('UnbHaar',x,Mfine,Mcoarse);
%        for Unbalanced Haar smoothing with integrals of scaling functions in
%        Mfine. size(Mfine) and size(Mcoarse) are optional, but if provided,
%        they must be equal to size(x)
%    H = constructprefilter('kernel',x,kerneltype,bandwidth,resolutionlevel);
%    H = constructprefilter('kernel',x,...
%          'type',kerneltype,'bandwidth',hJ,'level',j);
%    H = constructprefilter('wavelet',x,ht,sht)
%  See also
%    help RWT_2GprefiltPpolintUpdate


if length(x) == 1, 
   n = x; x = (0:n-1)/(n-1);
else
   n = length(x);
end

nvarargin = length(varargin);

switch H,
case {'none','I','eye','ident','identical'},
   % H = sparse(eye(n));
   H = 1;
case {'Haar','haar'},
   H = sparse(eye(n)+diag(ones(n-1,1),1))/2; H(n,n) = 1;
case {'Unbalanced Haar','UnbHaar'},
   Mfine = 0;
   if nvarargin > 0, if isreal(varargin{1}),
      Mfine = varargin{1};
   end, end
   if size(Mfine,1) < n | size(Mfine,2) < 1,
   % then varargin{1} was probably not meant as Mfine, so use default values
      Mfine = initmom(x);
   end
   Mfine = Mfine(1:n,1);
   if nvarargin > 1,
      Mcoarse = varargin{2};
   else
      Mcoarse = Mfine + [Mfine(2:n); 0];
      % averages with left evens and right odds. For averages with right evens
      % and left odds, use the following line:
      % Mcoarse = Mfine + [0; Mfine(1:n-1)];
   end
   if size(Mfine,1) ~= n, Mfine = Mfine'; end
   if size(Mcoarse,2) ~= size(Mfine,2), Mcoarse = Mcoarse'; end
   if size(Mcoarse,1) ~= size(Mfine,1), 
     error('Mcoarse must have the same size as Mfine')
   end
   % if size(Mcoarse,1) == ceil(n/2),
   %   Mcc = 
    % to be worked out...
   Mfine = Mfine(:,1); Mcoarse = Mcoarse(:,1);
   H = sparse(diag(Mfine./Mcoarse)+diag(Mfine(2:n)./Mcoarse(1:n-1),1));
   % the following line (in comment) for computing averages with right evens
   % and left odds:
   % H = sparse(diag(Mfine./Mcoarse)+diag(Mfine(1:n-1)./Mcoarse(2:n),-1));
   %
   % In the simple Unbalanced Haar transform, 
   % Mcoarse(k) = sum(Mfine(2*k-1:2*k)), but this 
   % need not be the case with general dual lifting steps.
   % In the simple Unbalanced Haar transform, Mfine can be computed from the
   % locations x at this scale, without knowing the finer scales. In general
   % lifting Mfine needs to be updated and transfered from scale to scale.
case {'kernel','Kernel'}
   k = 1;
   hJ = NaN;
   side = NaN;
   j = 0;
   while k<nvarargin,
      varname = varargin{k};
      if ischar(varname),
         switch varname,
         case {'bandwidth','width','scale'},
            hJ = varargin{k+1}; k = k+1;
         case {'type','kerneltype'},
            kerneltype = varargin{k+1}; k = k+1;
         case {'level','resolution'},
            j = varargin{k+1}; k = k+1;
         case {'side','Side'},
            side = varargin{k+1}; k = k+1;
         case {'left','right','L','R','both','twoside','two','center','Center'},
            side = varname;
         otherwise
            if exist('kerneltype') ~= 1,
               kerneltype = varname;
            end
         end
      elseif isnumeric(varargin{k}),
         if isnan(hJ), hJ = varargin{k};
         elseif j == 0, j = varargin{k};
         elseif isnan(side), side = varargin{k};
         end
      end
      k = k+1;
   end
   if k == nvarargin & isnumeric(varargin{k}),
      if isnan(hJ), hJ = varargin{k};
      elseif j == 0, j = varargin{k};
      elseif isnan(side), side = varargin{k};
      end
   end
   if isnan(side), side = 0; end
   N = length(x);
   hj = hJ*2^j;
   H = kernelestimation01(row(x),row(x),sparse(eye(N)),hj,kerneltype,side)';
   H = sparse(H);
case {'wavelet','wavfilt','2scale'}
   ht = [1 1]/2; sht = 0;
   if nvarargin > 0, ht = varargin{1}; end
   if nvarargin > 1, sht = varargin{2}; end
   % two lines below are correct but call of toeplitz is inefficient w.r.t.
   % memory use:
   % nht = length(ht); hht = zeros(1,n); hht(1:nht) = ht; hht = shift(hht,-sht);
   % H = toeplitz(hht,shift(hht(n:-1:1),-1));
   H = sparse(zeros(n,n));
   for k = 1:length(ht),
      o = 1-k-sht;
      H = H + sparse(diag(ones(1,n-abs(o))*ht(k),o));
      % next line for periodic extension only
      if o ~= 0,
         H = H + sparse(diag(ones(1,abs(o))*ht(k),-sign(o)*(n-abs(o))));
      end
   end
end

