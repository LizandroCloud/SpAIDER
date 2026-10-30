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
ts0 = 0.0;
%u = [ 25 23 21 18 24 23 22 20 19 28];  % serao variaveis de decisao do otimizador...
u = [ 0.009 0.0010 0.0050 0.0018 0.0015 0.0019 0.0016 0.00017 0.0019 0.00155];  % serao variaveis de decisao do otimizador...
ts = [ 111110.5 111110.3 110.4 1 0.2 0.7 0.5 0.2 0.12 0.8];  % serao valores calculados pelo otimizador...
hy = hyconnect;   %Conectar Hysys ao activeX controller.
hysys = hy.SimulationCases.Open([cd,'\MLBDEMO.usc']); %Abrindo a simulação
% hyissolving(hysys);
% hysolvertoggle(hysys);
% hyissolving(hysys);
h = get(hy.activeDocument.Flowsheet.Operations,'Item','SPRDSHT-1');  %
c = hycell(h, {'B1', 'B2', 'B3'});
times = hyvalue(c);
hyset(c{1},0);  %Change value of cell in Hysys spreadsheet.
hyintegtoggle(hy); %Toggle integrator on/off (dynamics mode).

for i=1:10

hysys.Activate  %Toggle solver on/off
hyissolving(hysys) %Check if solver is running.
 hyvalue(c)  %Get value of Hysys spreadsheet cell.
 hyset(c{2},ts(i)); %Change value of cell in Hysys spreadsheet.
 hyset(c{3},u(i)); %Change value of cell in Hysys spreadsheet.
 hyintegtoggle(hy); %Toggle integrator on/off (dynamics mode).


%hysolvertoggle(hysys);
%hyissolving(hysys);
%hyisintegrating(hysys) %Check if integrator is running (dynamics mode).
hyintegtime(hy)
hyvalue(c{3})
% hyrelease(hysys);
x(i) = hyvalue(c{3});
y(i) = hyvalue(c{2});
t(i) = ts(i) + ts0;
ts0 = t(i);
end
figure(1)
plot(t,x); hold on;
plot(t,x,'o');
figure(2)
plot(t,y); hold on;
plot(t,y,'o');
