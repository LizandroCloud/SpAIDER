
function y = IWT_2GPpolintU(x,w,nt,n,M);

%  FWT_2GPpolintU -- 1D, 2nd gen. inv. wav. transf.; polynomial interpol.
%                    prediction; prediction first
%   Usage
%     y       = FWT_2GPpolintU(x,w,nt,n,M);
%   Inputs
%     x      grid-locations of observations (time, space)
%     w      wavelet coefficients
%     nt     dual vanishing moments -- polynomial prediction order
%     n      primal vanishing moments -- number of taps in update filter
%     M      number of resulution levels
%   Outputs
%     y      vector of size(w) with reconstructed signal/function/observations
%   Description
%     Applies repeated lifting filter bank with one prediction step and one
%     update step. The prediction is based on polynomial interpolation, the
%     update uses n taps (i.e., n nonzeros, i.e., n degrees of freedom)
%   Note
%   See also
%     help FWT_2GPpolintU


transposex = 0;
if size(x,2) > size(x,1),
   x = x';
end
transposew = 0;
if size(w,1) ~= size(x,1),
   w = w';
end

moments = cell(M);
grids = cell(M);
if n > 0,
   Mfine = initmom(x,(0:n-1));
else
   Mfine = zeros(length(x),0);
end
xfine = x;

for j = 1:M
   moments{j} = Mfine;
   grids{j} = xfine;
   [even, odd] = splitevenodd(xfine);
   xeven = xfine(even);
   xodd = xfine(odd);
   Ne = length(even);
   Ie = sparse(eye(Ne));
   P = sparse(polynomialprediction(xeven,Ie,xodd,nt));
   Meven = Mfine(even,1:n); Modd = Mfine(odd,1:n);
   Mcoarse = Meven + P'*Modd;
   Mfine = Mcoarse;
   xfine = xeven;
end

y = w;
yeven = y(1:Ne,:);
for j = M:-1:1,
   xfine = grids{j};
   [even, odd] = splitevenodd(xfine);
   xeven = xfine(even);
   xodd = xfine(odd);
   Mfine = moments{j};
   Nj = length(xfine);
   Ne = length(xeven);
   No = Nj-Ne;

   yodd = w(Ne+1:Nj,:);
   Ie = sparse(eye(Ne));
   P = sparse(polynomialprediction(xeven,Ie,xodd,nt));
   % UPDATE
   Modd = Mfine(odd,:);
   if n > 0,
      yeven = yeven - updateprimalmoments(xodd,yodd,xeven,Modd,Mcoarse);
   end
   Mcoarse = Mfine;
   % PREDICTION %
   yodd = yodd + P*yeven;
   % MERGE %
   y(even,:) = yeven;
   y(odd,:) = yodd;
   yeven = y(1:Nj,:);
end

if transposew, y = y'; end
