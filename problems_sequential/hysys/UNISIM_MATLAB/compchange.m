% % Composition change in a stream of HYSYS using ActiveX server 
% % This is done by making a macro in HYSYS and playing it in MATLAB
% % by: Maryam Ghadrdan (maryam.ghadrdan@gmail.com) & Hamid Mehdizadeh
% % (hamid.mehdizadeh@gmail.com)
% % 27-Feb-2011
% 
% clc
% clear all
% fid = fopen('changecomp1.scp', 'r+');
% fseek(fid, 188, 'bof');
% fwrite(fid, num2str(7));
% fclose(fid);
% type('changecomp1.scp')

% serv = actxserver('hysys.application.2006');
% hyCase = serv.ActiveDocument;
% % hyCase.Application.PlayScript([pwd],'\changecomp1.SCP')
% Feed = hyCase.Flowsheet.MaterialStreams.Item('Feed');
% hyCase.Application.PlayScript('C:\Users\lades\Documents\MATLAB\changecomp1.scp')
for i=1:10
hysys = hyconnect;
sprd = hyspread(hysys, 'SPRDSHT-1');
c = hycell(sprd, {'B1', 'B2', 'B3', 'B4'});
hysolvertoggle(hysys);
hyissolving(hysys);
 hyvalue(c)
 hyset(c{3},22+i);
hysolvertoggle(hysys);
hyintegtoggle(hysys);
hyissolving(hysys);
hyisintegrating(hysys);
hyintegtime(hysys);
hyvalue(c{2})
% hyrelease(hysys);
x(i) = hyvalue(c{4});
end
plot(x,'o');