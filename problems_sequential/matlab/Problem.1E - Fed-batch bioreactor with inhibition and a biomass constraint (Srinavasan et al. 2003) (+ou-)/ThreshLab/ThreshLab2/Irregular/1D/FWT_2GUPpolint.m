
function [w, Ne] = FWT_2GUPpolint(x,y,nt,n,M);

%  FWT_2GPpolintU -- 1D, 2nd gen. frwd. wav. transf.; polynomial interpol.
%                    prediction; update first
%   Usage
%     [w, Ne] = FWT_2GUPpolint(x,y,nt,n,M);
%     w       = FWT_2GUPpolint(x,y,nt,n,M);
%   Inputs
%     x      grid-locations of observations (time, space)
%     y      observed values
%     nt     dual vanishing moments -- polynomial prediction order
%     n      primal vanishing moments -- number of taps in update filter
%     M      number of resulution levels
%   Outputs
%     w      vector of size(y) with wavelet coefficients
%     Ne     (optional) number of scaling coefficients (dependent on M)
%   Description
%     Applies repeated lifting filter bank with one prediction step and one
%     update step. The prediction is based on polynomial interpolation, the
%     update uses n taps (i.e., n nonzeros, i.e., n degrees of freedom)
%   Note
%     Between update and prediction, the updated evens are rescaled
%   See also
%     help IWT_2GUPpolint


transposex = 0;
if size(x,2) > size(x,1),
   x = x';
end
transposey = 0;
if size(y,1) ~= size(x,1),
   y = y';
end

if n > 0,
   Mfine = initmom(x,(0:n-1));
end

xeven = x;
yeven = y;
w = zeros(size(y));
for j = 1:M

%%%%%%%  scale = j

  % SPLIT %
  [even, odd] = splitevenodd(x);
  yodd  = y(odd,:);
  yeven = y(even,:);
  xodd = x(odd);
  xeven = x(even);
  Ne = length(even); No = length(odd);
  if n>0,
     Meven = Mfine(even,1:n); Modd = Mfine(odd,1:n); 
  end
  % UPDATE
  if n > 0,
     yeven = yeven + updateprimalmoments(xodd,yodd,xeven,Modd,Meven);
     e1 = ones(size(xeven)); o1 = ones(size(xodd));
     rescale = e1 + updateprimalmoments(xodd,o1,xeven,Modd,Meven);
     D = sparse(diag(rescale));
     yeven = D\yeven;
  end
  % PREDICTION MATRIX %
  I = sparse(eye(Ne));
  P = sparse(polynomialprediction(xeven,I,xodd,nt));
  % PREDICTION %
  yodd = yodd - P*yeven;
  % MOMENTS %
  if n > 0,
     Mfine = D*Meven;
  % in case U = 0, moments are different, but they aren't used, so computation
  % is not necessary:
  % else
  %     Mfine = Meven + P'*Modd;
  end
  % FILLING IN %
  x = xeven;
  y = yeven;
  w(Ne+1:Ne+No,:) = yodd;

end
w(1:Ne,:) = yeven;

if transposey, w = w'; end
