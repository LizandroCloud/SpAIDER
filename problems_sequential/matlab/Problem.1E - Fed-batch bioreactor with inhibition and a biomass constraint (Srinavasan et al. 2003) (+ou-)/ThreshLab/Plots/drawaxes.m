
function void = drawaxes(varargin);

% drawaxes - ThreshLab/Plots -- 
%   Usage
%     drawaxes
%   Inputs
%     void
%   Outputs
%     void
%   Description
%     Draws arrowed axes on current figure


arrowlength = 0.07;
arrowwidth = 0.02;
holdon = ishold;
axy = axis;
x0 = axy(1); x1 = axy(2); y0 = axy(3); y1 = axy(4);
if x0*x1 > 0, Ypos = x0; else, Ypos = 0; end
if y0*y1 > 0, Xpos = y0; else, Xpos = 0; end

Dx = x1-x0; Dy = y1-y0;
axy(1) = min(axy(1),Ypos-Dx*arrowwidth);
axy(2) = max(axy(2),Ypos+Dx*arrowwidth);
axy(3) = min(axy(3),Xpos-Dy*arrowwidth);
axy(4) = max(axy(4),Xpos+Dy*arrowwidth);



hold on
plot(axy(1:2)-[0 1]*Dx*arrowlength,Xpos*[1 1],'linewidth',3,'color','k')
plot(Ypos*[1 1],axy(3:4)-[0 1]*Dy*arrowlength,'linewidth',3,'color','k')

Dx = axy(2)-axy(1);
Dy = axy(4)-axy(3);
fill(axy(2)-Dx*arrowlength*[1 0 1],Xpos+arrowwidth*Dy*[-1 0 1],'k')
fill(Ypos+arrowwidth*Dx*[-1 0 1],axy(4)-Dy*arrowlength*[1 0 1],'k')


axis(axy)
if ~holdon, hold off, end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
