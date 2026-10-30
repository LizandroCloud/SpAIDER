function [tspan,y] = model_unisim_o(ti,dt,u,ts,ks,hy,hysys)

% Software para Comunicação UNSIM x MATLAB
% NEO-UFF - Universidade Federal Fluminense
% Última modificação - 17 de Novembro de 2017 - 16 h 24 min

% u=86;
% Comandos usados para ler e escrever informações no UNiSIM. 
% u = [ 0.009 0.0010 0.0050 0.0018 0.0015 0.0019 0.0016 0.00017 0.0019 0.00155];  % serao variaveis de decisao do otimizador...
% ts = [ 0.5 0.3 0.4 1 0.2 0.7 0.5 0.2 0.12 0.8];  % serao valores 
% hy = hyconnect;   %Conectar Hysys ao activeX controller.
% hysys = hy.SimulationCases.Open([cd,'\gstrt.usc']); %Abrindo a simulação
% hysys.Visible = true; %Visibilidade
b = get(hy.activeDocument); %importanto dados do objeto principal
c = get(b.Flowsheet); %importanto dados do flowsheet
d = get(c.Operations); %importanto dados das operações
% d.Names; %importanto os nomes dos equipamentos
h = get(hy.activeDocument.Flowsheet.Operations,'Item','TABLE');  %
g = get(hy.activeDocument.Flowsheet.Operations,'Item','TAB');  %
i = get(h); %importanto o equipamento V-100

A = hycell(h, {'A2', 'A3','A3', 'A4','A5', 'A6','A7'});
Ao = hyvalue(A);
C = hycell(h, {'C4'});
TAB = hycell(g, {'B1', 'B2'});
a1 = hyvalue(C);  %Get value of Hysys spreadsheet cell.
times = hyvalue(TAB);
% volume = i.VesselVolume.value; %lendo o valor do volume do vaso
% i.VesselVolume.value = 5; %escrevendo o valor do volume do vaso


if ks==1
    hyset(TAB{1},ts(ks));  %Change value of cell in Hysys spreadsheet.
    tt = linspace(ts(ks),ti,dt);
    to=0;
else
    hyset(TAB{1},ts(ks));
    tt = linspace(ts(ks),ti,dt);
    to=ts(ks);
end
t_final = (ti-to)/10;
tx=ts(ks);

hyset(TAB{1},0);  %Change value of cell in Hysys spreadsheet.
hyset(TAB{2},t_final);

for i=1:dt
    
    hyset(C{1},u);
%     hyset(TAB{1},0);  %Change value of cell in Hysys spreadsheet.
%     hyset(TAB{2},t_final);
    hysolvertoggle(hysys) ; %Toggle solver on/off
    hyintegtoggle(hysys);
    hysolvertoggle(hysys);
    resp(i,1:7) = hyvalue(A);
    tspan(i) = tx+t_final;
    tx = tspan(i);
end
tspan = [to tspan];
y= [Ao; resp];
end