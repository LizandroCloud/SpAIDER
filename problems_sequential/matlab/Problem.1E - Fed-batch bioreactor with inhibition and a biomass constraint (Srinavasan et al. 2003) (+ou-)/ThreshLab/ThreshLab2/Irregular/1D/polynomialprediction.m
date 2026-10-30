
function P = polynomialprediction(xeven,yeven,xodd,nt)

% polynomialprediction - ThreshLab2/Irregular/1D/ --- polynomial prediction
% lifting step
%   Usage
%      yodd = polynomialprediction(xeven,yeven,xodd,nt);
%      P    = polynomialprediction(xeven,I,xodd,nt);
%   Inputs
%      xeven  x-coordinates of even indexed locations
%      yeven  y-values in even indexed locations
%      I      matrix (e.g., identity matrix of size = length(xeven) )
%      xodd   x-coordinates of odd indexed locations
%      nt     number of dual vanishing moments (= order of polynomial
%             prediction, i.e.: 2 = linear; 4 = cubic, etc.)
%   Outputs
%      yodd   predicted values in odd indexed locations
%      P      for I = eye(length(xeven)), P is the prediction operation matrix
%   Description
%      This routine evaluates an interpolating polynomial in each point of xodd.
%      The interpolating polynomial may be different for every point in xodd.
%      if nt == 2, this routine coincides with polyline except for points
%                  outside [min(xeven),max(xeven)]
%   See Also
%      help interpolationneville
%      help polyline

Ne = length(xeven);
Noo = length(xodd);
[xodd oo] = sort(xodd);
% we assume xeven and xodd are row vectors. If xodd turns out to be column
% vector, then a final tranpose is applied to P
Ptranspose = 0;
if size(xeven,1) > 1, xeven = xeven'; Ptranspose = 1; end
if size(xodd,1) > 1, xodd = xodd'; Ptranspose = 1; end
[r s] = size(yeven);
if s ~= Ne, 
   yeven = yeven';
   [r s] = size(yeven);
   if s ~= Ne, error('dimension mismatch inputs xeven & yeven'), end
end

P = zeros(r,Noo);
if (Ne < nt) 

  nt = Ne;
  P = interpolationneville(xodd,xeven,yeven);

else

  lefteven = rankoddsbetweenevens(xeven,xodd);
  No = length(xodd); % maximum value is Ne+1

  nl = ceil(nt/2); nr = nt-nl;
  % if nt is odd, then the prediction will be based on nl even points on the
  % left and nl-1 on the right. As the first point in a split is always even,
  % there is a slight tendency to have more points on the left available
  for e = 1:nl-1, o = (lefteven==e); if any(o),
    P(:,o) = interpolationneville(xodd(o),xeven(1:nt),yeven(:,1:nt));
  end, end
  for e = nl:Ne-nr, o = (lefteven==e); if any(o),
    P(:,o) = ...
    interpolationneville(xodd(o),xeven(e-nl+1:e+nr),yeven(:,e-nl+1:e+nr));
  end, end
  for e = Ne-nr+1:No, o = (lefteven==e); if any(o),
    P(:,o) = ...
    interpolationneville(xodd(o),xeven(Ne-nt+1:Ne),yeven(:,Ne-nt+1:Ne));
  end, end

end
if Ptranspose, P = P'; end
