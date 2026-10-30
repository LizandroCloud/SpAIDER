
function ii = rankoddsbetweenevens(xeven,xodd);

% rankoddsbetweenevens - ThreshLab2/Irregular/1D/ --- index manipulation
%   Usage
%      ii = rankoddsbetweenevens(xeven,xodd);
%   Inputs
%      xeven  x-coordinates of even indexed locations, in ascending order
%      xodd   x-coordinates of odd indexed locations, in arbitrary order
%   Outputs
%      ii     integer vector of length=length(xodds) with index into xeven
%             xeven(ii(j)) is largest "even" point smaller than or equal to
%             xodd(j)
%             ii(j) = 0 if no point in xeven is below xodd(j)
%   Note
%      This function is used in polynomialprediction
%   See also
%      help polynomialprediction

  No = length(xodd);
  Ne = length(xeven);
  ii = Ne*ones(1,No);
  for i = 1:Ne,
     ii(xodd < xeven(i)) = ii(xodd < xeven(i))-1;
  end

