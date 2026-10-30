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

% serv = actxserver('hysys.application');
% hyCase = serv.ActiveDocument;
% % hyCase.Application.PlayScript([pwd],'\changecomp1.SCP')
% Feed = hyCase.Flowsheet.MaterialStreams.Item('Feed');
% hyCase.Application.PlayScript('C:\Users\lades\Documents\MATLAB\changecomp1.scp')
ts0 = 0.0;
%u = [ 25 23 21 18 24 23 22 20 19 28];  % serao variaveis de decisao do otimizador...
u = [ 0.009 0.0010 0.0050 0.0018 0.0015 0.0019 0.0016 0.00017 0.0019 0.00155];  % serao variaveis de decisao do otimizador...
ts = [ 0.5 0.3 0.4 1 0.2 0.7 0.5 0.2 0.12 0.8];  % serao valores calculados pelo otimizador...
hy = hyconnect;   %Connect to Hysys application as a activeX controller.
hy.Visible = true;
hysys = hy.SimulationCases.Open([cd,'\eb2dyn.usc']);
hysys.Visible = true;
b = get(hy.activeDocument);
c = get(b.Flowsheet);
d = get(c.Operations);
d.Names;
h = get(hy.activeDocument.Flowsheet.Operations,'Item','V-100');
i = get(h);
temp = i.VesselVolume.value;
% sprd = hyspread(hysys, 'SPRDSHT-1');  %Connect to Hysys spreadsheet.
% c = hycell(sprd, {'B1', 'B2', 'B3','B4'});
% hyset(c{1},0);  %Change value of cell in Hysys spreadsheet.
% hyintegtoggle(hy); %Toggle integrator on/off (dynamics mode).
for j=1:10

hysolvertoggle(hy)  %Toggle solver on/off
%hyissolving(hysys) %Check if solver is running.
%  hyvalue(c)  %Get value of Hysys spreadsheet cell.
%  hyset(c{2},ts(i)); %Change value of cell in Hysys spreadsheet.
%  hyset(c{3},u(i)); %Change value of cell in Hysys spreadsheet.
  %Toggle integrator on/off (dynamics mode).
% 
b = get(hy.activeDocument);
c = get(b.Flowsheet);
d = get(c.Operations);
d.Names;
h = get(hy.activeDocument.Flowsheet.Operations,'Item','V-100');
i = get(h);
% temp = i.VesselVolume.value;
%hysolvertoggle(hysys);
%hyissolving(hysys);
%hyisintegrating(hysys) %Check if integrator is running (dynamics mode).
hyintegtime(hy)
% hyvalue(c{3})
% hyrelease(hysys);
% x(i) = hyvalue(c{3});
% y(i) = hyvalue(c{4});
hyintegtoggle(hy);
t(j) = ts(j) + ts0;
ts0 = t(j);
temp(j) = i.VesselTemperatureValue;
end
figure(1)
plot(t,temp); hold on;
% plot(t,x,'o');
% figure(2)
% plot(t,y); hold on;
% plot(t,y,'o');
