clear all
close all
% clc

global k alfa epsilon kla T R P x_o2 V_g V_L F_ar c_m c_g0 t_i t_f

t_i = 8;
t_f = 8.01;

load('iJO1366_PM_asp_plasm_novo_1.mat');
modelo=iJO1366_PM_asp_plasm_novo_1;

modelo = changeRxnBounds(modelo,'EX_glc(e)',0,'b');
modelo = changeRxnBounds(modelo,'EX_glyc(e)',(-7),'l'); 
modelo = changeRxnBounds(modelo,'EX_glyc(e)',0,'u'); 
modelo = changeRxnBounds(modelo,'NADTRHD',4,'u'); % NAD transhydrogenase
modelo = changeRxnBounds(modelo,'G6PDH2r',3.5,'u'); % glucose 6 phosphate 
% dehydrogenase - via das pentoses

k=0.8;
alfa=50;
epsilon=1e-6;
kla=41.3;
T=37;
R=0.082; %atm L/mol K
P=1; %atm
x_o2=0.21;
V_g=1;%L
V_L=1.5;%L
F_ar=90;%L/h

%Fator de conversão enzimático: 300
c=300;

ID_o2 = findMetIDs(modelo,'o2[e]'); 

% Conc. de o2 na vazão de alimentação de ar (c_g0)
c_g0=P*modelo.PM(ID_o2)*x_o2/(R*(T+273.15));

% Conc. de saturação de o2 na fase líquida
c_m=3.62e-06*T^2-3.312e-04*T+0.01466217;

%Condições iniciais

V0=1.5;
X0=0.05;
G0=30;
Ac0=0;
Lac0=0;
L0=0;
Asp0=0;
Plasm0=0;
C_L0=0.0072;
C_G0=c_g0;

x=[V0 X0 G0 Ac0 Lac0 L0 Plasm0 Asp0 C_L0 C_G0];

% [t,x]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo),[0 25],x0);

S=0.1:0.1:1;
m=length(S);
x(3)=0.4;
x(9)=0;
for i=1:m
 x(4)=S(i);
 v=pseudo_steady_state_v(0,x,modelo);
 [ID,f]=fluxos(modelo,v);
 asp(i)=f(5);
 mi(i)=v(8);
end

figure(1)
plot(S,mi')
xlabel('Substrato')
ylabel('mi')

figure(2)
plot(S,asp')
xlabel('Substrato')
ylabel('asp')

