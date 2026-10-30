clear all
close all
clc

load('iJO1366_PM_asp_plasm_novo_1.mat')
 
modelo=iJO1366_PM_asp_plasm_novo_1;

v_plasm_max=2e-6;
% v_plasm_max=1.5767E-4;
% v_plasm_max=1E-5;
epsilon=linspace(0,v_plasm_max,50);

for i=1:length(epsilon)
[v,obj,erro]=fba2_multi(modelo,0.8,0.3899,50,epsilon(i),2);
v_plasm(i)=v(2587);
v_asp(i)=v(2585);
fun_obj(i)=obj;
v_o2(i)=v(252);
v_acetato(i)=v(512);
v_lactate(i)=v(1622);
end

subplot(2,3,1),plot(epsilon,v_asp)
xlabel('epsilon')
ylabel('Fluxo da asparaginase (mmol/g h)')
subplot(2,3,2),plot(epsilon,fun_obj)
xlabel('epsilon')
ylabel('Função objetivo')
subplot(2,3,3),plot(epsilon,v_acetato)
xlabel('epsilon')
ylabel('Fluxo do acetato (mmol/g h)')
subplot(2,3,4),plot(epsilon,v_lactate)
xlabel('epsilon')
ylabel('Fluxo do lactato (mmol/g h)')
ylim([-1 1])
subplot(2,3,5:6),plot(epsilon,v_o2)
xlabel('epsilon')
ylabel('Fluxo do o2 (mmol/g h)')

%valor máximo experimental do fluxo de plasmídeo: 1.27e-6
%2009_Williams et al_Biotechnology and Bioengineering