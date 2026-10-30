% Software para Comunicação UNSIM x MATLAB
% NEO-UFF - Universidade Federal Fluminense
% Última modificação - 17 de Novembro de 2017 - 16 h 24 min


ts0 = 0.0; % tempo inicial
ts = [ 0.5 0.3 0.4 1 0.2 0.7 0.5 0.2 0.12 0.8];  % dts para integração
hy = hyconnect;   %Conectar Hysys ao activeX controller.
hysys = hy.SimulationCases.Open([cd,'\dyntut3.hsc']); %Abrindo a simulação
hysys.Visible = true; %Visibilidade
sprd = hysys.Flowsheet.Operations.Item('SPRDSHT-1');
a = sprd.Cell('B2').CellValue
% sprd.Cell('B2').CellValue = 0.75
% hyintegtoggle(hy);
c = hycell(sprd, {'B1', 'B2', 'B3', 'B4'});
hysys.Visible = true; %Visibilidade
b = get(hy.activeDocument); %importanto dados do objeto principal
c = get(b.Flowsheet); %importanto dados do flowsheet
d = get(c.Operations); %importanto dados das operações
% d.Names; %importanto os nomes dos equipamentos
% h = get(hy.activeDocument.Flowsheet.Operations,'Item','V-100');  %importanto o equipamento V-100
% i = get(h); %importanto o equipamento V-100
% volume = i.VesselVolume.value; %lendo o valor do volume do vaso
% i.VesselVolume.value = 5; %escrevendo o valor do volume do vaso
% b = get(hy.activeDocument);
% c = get(b.Flowsheet);
% d = get(c.Operations);
% d.Names;
% hysolvertoggle(hy)  %Toggle solver on/off
hyintegtime(hy)
hyintegtoggle(hy);
sprd.Cell('B1').CellValue = 0;
for j=1:10 %testando uma iteração

sprd.Cell('B2').CellValue = ts(j);
sprd.Cell('B3').CellValue = ts(j);
% b = get(hy.activeDocument);
% c = get(b.Flowsheet);
% d = get(c.Operations);
% d.Names;
% h = get(hy.activeDocument.Flowsheet.Operations,'Item','V-100');
% i = get(h);
hyintegtime(hy)
hyintegtoggle(hy);
t(j) = ts(j) + ts0;
ts0 = t(j);
temp(j) = i.VesselTemperatureValue;
end
%% 
plot(t, temp)
