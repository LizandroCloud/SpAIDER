
function void = verticalboxplot(X,varargin)

% /ThreshLab/General/verticalboxplot -- produces vertical (i.e., standing) 
%                                       boxplots
%  Usage
%    verticalboxplot(X,'shaded','position',p,'width',w,plotoptions)
%  Inputs
%    X               vector, matrix, or cell array of data
%    'position', p   vector of positions/locations of the boxplots (optional)
%    'width',w       vector of widths of the boxplots (optional)
%    'shaded'        makes central box shaded (optional) 
%  Outputs
%    <none>
%  Description
%     This routine plots a single box if X is a vector, or multiple boxes if X
%     is a cell array (one box for each cell) or a matrix (one box for each
%     column of X).
%     The boxes are plotted at the locations indicated in vector locations
%     if no locations are given, then boxes are plotted at equidistant
%     locations
%   Note
%   See also
%     help horizontalboxplot

po = 1;
ia = 1;
fillboxes = 0;
while ia < nargin,
  switch varargin{ia}
  case {'shaded','shad'},
     fillboxes = 1; ia = ia+1;
  case {'position','pos'}
     pos = varargin{ia+1}; ia = ia+2;
  case {'width'}
     width = varargin{ia+1}; ia = ia+2;
  otherwise
     plotoptions{po} = varargin{ia}; ia = ia+1; po = po+1;
  end
end
if po==1, plotoptions = cell(0,0); end

hh = ishold;

if ~iscell(X),
   rc = size(X);
   if min(rc) == 1, X = column(X); end
   cc = size(X,2);
   for c = 1:cc,
      x = X(:,c);
      XX{c} = x;
   end
   X = XX;
end
if exist('pos') ~= 1,
   pos = (1:length(X));
end
if exist('width') ~= 1,
   width = ones(1,length(X))/5;
end
if length(pos) ~= length(X),
   pos = (1:length(X));
end
if length(width) ~= length(X),
   width = ones(1,length(X))/5;
end
for c = 1:length(X),
   x = X{c};
   h = width(c);
   y = pos(c);

   x2 = median(x);
   xx = x(x<x2);
   x1 = median(xx);
   xx = x(x>x2);
   x3 = median(xx);
   IQR = x3 - x1;
   xbar = mean(x);
   w2 = min(max(x),x3+IQR*1.5);
   w1 = max(min(x),x1-IQR*1.5);
   outliers = [x(x<w1); x(x>w2)];
   plot([y-h y+h y+h y-h y-h],[x1 x1 x3 x3 x1],plotoptions{:})
   hold on
   if fillboxes,
      fillcolor = get(findobj(gcf,'type','line'),'color');
      if iscell(fillcolor),
         fillcolor = fillcolor{1};
      end
      fillcolor = 1-(1-fillcolor)/4;
      fill([y-h y+h y+h y-h y-h],[x1 x1 x3 x3 x1],fillcolor)
      plot([y-h y+h y+h y-h y-h],[x1 x1 x3 x3 x1],plotoptions{:})
   end
   plot([y-h y+h],[x2 x2],plotoptions{:})
   plot([y y],[x3 w2],plotoptions{:})
   plot([y-h/2 y+h/2],[w2 w2],plotoptions{:})
   plot([y y],[w1 x1],plotoptions{:})
   plot([y-h/2 y+h/2],[w1 w1],plotoptions{:})
   plot(y*ones(size(outliers)),outliers,'+',plotoptions{:})
   plot(y,xbar,'*',plotoptions{:})
end

if ~hh, hold off, end
