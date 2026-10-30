function [tspan,y] = model_hysys_o(ti,dt,u,ts,ks,hy,hysys)

ts0 = 0;
hy = hyconnect;   %Conectar Hysys ao activeX controller.
hysys = hy.SimulationCases.Open([cd,'\dyntut3.hsc']); %Abrindo a simulação
hysys.Visible = true; %Visibilidade
sprd = hysys.Flowsheet.Operations.Item('SPRDSHT-1');
hysys.Visible = true; %Visibilidade
hyintegtime(hy)
hyintegtoggle(hy);
sprd.Cell('B1').CellValue = 0;
for j=1:length(ts)-1 %testando uma iteração
sprd.Cell('C5').CellValue = u(j) + rand(1)*10;
% sprd.Cell('B2').CellValue = ts(j);
sprd.Cell('B3').CellValue = ts(j);
hyintegtime(hy)
hyintegtoggle(hy);
t(j+1) = ts(j+1) + ts0;
ts0 = t(j+1);
x1(j) = sprd.Cell('D1').CellValue;
x2(j) = sprd.Cell('D2').CellValue;

end



end