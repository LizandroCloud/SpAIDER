% Software para Comunicação UNSIM x MATLAB
% NEO-UFF - Universidade Federal Fluminense
% Última modificação - 17 de Novembro de 2017 - 16 h 24 min


% Comandos usados para ler e escrever informações no UNiSIM. 
u = [ 0.009 0.0010 0.0050 0.0018 0.0015 0.0019 0.0016 0.00017 0.0019 0.00155];  % serao variaveis de decisao do otimizador...
ts = [ 0.5 0.3 0.4 1 0.2 0.7 0.5 0.2 0.12 0.8];  % serao valores 
hy = hyconnect;   %Conectar Hysys ao activeX controller.
hysys = hy.SimulationCases.Open([cd,'\gstrt.usc']); %Abrindo a simulação
hysys.Visible = true; %Visibilidade
b = get(hy.activeDocument); %importanto dados do objeto principal
c = get(b.Flowsheet); %importanto dados do flowsheet
d = get(c.Operations); %importanto dados das operações
% d.Names; %importanto os nomes dos equipamentos
h = get(hy.activeDocument.Flowsheet.Operations,'Item','TABLE');  %
g = get(hy.activeDocument.Flowsheet.Operations,'Item','TAB');  %
i = get(h); %importanto o equipamento V-100

A = hycell(h, {'A2', 'A3','A3', 'A4','A5', 'A6','A7'});
C = hycell(h, {'I3'});
TAB = hycell(g, {'B1', 'B2'});
a1 = hyvalue(C)  %Get value of Hysys spreadsheet cell.
times = hyvalue(TAB);
% volume = i.VesselVolume.value; %lendo o valor do volume do vaso
% i.VesselVolume.value = 5; %escrevendo o valor do volume do vaso
hyset(C{1},0.13);
hyset(TAB{1},0);  %Change value of cell in Hysys spreadsheet.
hyset(TAB{2},10);
hysolvertoggle(hysys)  %Toggle solver on/off
hyintegtoggle(hysys)
hysolvertoggle(hysys)
a=2;
