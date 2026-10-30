%clear all
close all
% clc

global k alfa epsilon kla T R P x_o2 V_g V_L F_ar c_m c_g0 t_i t_f

t_i = 8;
t_f = 8.01;

load('iJO1366_PM_asp_plasm_novo_2.mat');
modelo=iJO1366_PM_asp_plasm_novo_2;

modelo = changeRxnBounds(modelo,'EX_glc(e)',0,'b');
modelo = changeRxnBounds(modelo,'EX_glyc(e)',(-7),'l'); 
modelo = changeRxnBounds(modelo,'EX_glyc(e)',0,'u'); 
modelo = changeRxnBounds(modelo,'NADTRHD',4,'u'); % NAD transhydrogenase
modelo = changeRxnBounds(modelo,'G6PDH2r',3.5,'u'); % glucose 6 phosphate 
% dehydrogenase - via das pentoses

k=0.8;
alfa=50;
epsilon=1e-6;
kla=41.3/10;
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

x0=[V0 X0 G0 Ac0 Lac0 L0 Plasm0 Asp0 C_L0 C_G0];

% [t,x]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo),[0 25],x0);

op=odeset('MaxStep',0.01);
[t,x]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo),[0:1:t_i],x0);
[ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo),[t_i:0.001:t_f],x(end,:),op);
x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
[ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo),[t_f:0.0001:11],x(end,:));
x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
[ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo),[11:0.00001:20],x(end,:));
x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
% [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo),[20:0.0001:24],x(end,:));
% x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];

p_o2=100*x(:,9)/c_m;
Enz=c*x(:,8); % [U/mL]

if 0
m=length(t);
for i=1:m
 v=pseudo_steady_state_v(t(i),x(i,:),modelo);
 mi(i)=v(8);
end

figure(1)
plot(t,mi')
xlabel('Tempo')
ylabel('mi')

figure(2)
plot(x(:,3),mi')
xlabel('Substrato')
ylabel('mi')
end

% 
% rxns=findRxnIDs(iJO1366_PM_asp_plasm,{'RPI','TKT1','TKT2','TALA','GLYCK',...
% 'PGM','ENO','PPCK','PPC','PDH','MDH','CS','AKGDH','ICDHyr','ACKr','PTAr'});
% for i=1:m
% mi(i)=v(i,8);
% RPI(i)=v(i,rxns(1));
% TKT1(i)=v(i,rxns(2));
% TKT2(i)=v(i,rxns(3));
% TALA(i)=v(i,rxns(4));
% GLYCK(i)=v(i,rxns(5));
% PGM(i)=v(i,rxns(6));
% ENO(i)=v(i,rxns(7));
% PPCK(i)=v(i,rxns(8));
% PPC(i)=v(i,rxns(9));
% PDH(i)=v(i,rxns(10));
% MDH(i)=v(i,rxns(11));
% CS(i)=v(i,rxns(12));
% AKGDH(i)=v(i,rxns(13));
% ICDHyr(i)=v(i,rxns(14));
% ACKr(i)=v(i,rxns(15));
% PTAr(i)=v(i,rxns(16));
% end
% 
% figure(9)
% plot(t,mi')
% xlabel('Tempo')
% ylabel('mi')
% figure(10)
% plot(t,RPI')
% xlabel('Tempo')
% ylabel('RPI')
% figure(11)
% plot(t,TKT1')
% xlabel('Tempo')
% ylabel('TKT1')
% figure(12)
% plot(t,TKT2')
% xlabel('Tempo')
% ylabel('TKT2')
% figure(13)
% plot(t,TALA')
% xlabel('Tempo')
% ylabel('TALA')
% figure(14)
% plot(t,GLYCK')
% xlabel('Tempo')
% ylabel('GLYCK')
% figure(15)
% plot(t,PGM')
% xlabel('Tempo')
% ylabel('PGM')
% figure(16)
% plot(t,ENO')
% xlabel('Tempo')
% ylabel('ENO')
% figure(17)
% plot(t,PPCK')
% xlabel('Tempo')
% ylabel('PPCK')
% figure(18)
% plot(t,PPC')
% xlabel('Tempo')
% ylabel('PPC')
% figure(19)
% plot(t,PDH')
% xlabel('Tempo')
% ylabel('PDH')
% figure(20)
% plot(t,MDH')
% xlabel('Tempo')
% ylabel('MDH')
% figure(21)
% plot(t,CS')
% xlabel('Tempo')
% ylabel('CS')
% figure(22)
% plot(t,AKGDH')
% xlabel('Tempo')
% ylabel('AKGDH')
% figure(23)
% plot(t,ICDHyr')
% xlabel('Tempo')
% ylabel('ICDHyr')
% figure(24)
% plot(t,ACKr')
% xlabel('Tempo')
% ylabel('ACKr')
% figure(25)
% plot(t,PTAr')
% xlabel('Tempo')
% ylabel('PTAr')
% aa=findRxnIDs(modelo,{'GLYTRS','ALATRS','VALTRS','LEUTRS','ILETRS',
% 'PHETRS','PROTRS','SERTRS','THRTRS','TYRTRS','CYSTRS','METTRS','LYSTRS',
% 'ARGTRS','HISTRS','TRPTRS','ASPTRS','GLUTRS','ASNTRS','GLNTRS'});
% 
% for i=1:m
% mi(i)=v(i,8);
% GLYTRS(i)=v(i,aa(1));
% ALATRS(i)=v(i,aa(2));
% VALTRS(i)=v(i,aa(3));
% LEUTRS(i)=v(i,aa(4));
% ILETRS(i)=v(i,aa(5));
% PHETRS(i)=v(i,aa(6));
% PROTRS(i)=v(i,aa(7));
% SERTRS(i)=v(i,aa(8));
% THRTRS(i)=v(i,aa(9));
% TYRTRS(i)=v(i,aa(10));
% CYSTRS(i)=v(i,aa(11));
% METTRS(i)=v(i,aa(12));
% LYSTRS(i)=v(i,aa(13));
% ARGTRS(i)=v(i,aa(14));
% HISTRS(i)=v(i,aa(15));
% TRPTRS(i)=v(i,aa(16));
% ASPTRS(i)=v(i,aa(17));
% GLUTRS(i)=v(i,aa(18));
% ASNTRS(i)=v(i,aa(19));
% GLNTRS(i)=v(i,aa(20));
% end

% [t,x]=ode45('pseudo_steady_state',[0 6],x0);

% legenda={'Biomassa(g/L)','Glicerol(g/L)','Ácido acético(g/L)'};
% plotStyle = {'b-','k:','r.','yo'};
% for i=2:4
%     plot(t,x(:,i),plotStyle{i})
%     xlabel('Tempo')
%     legendInfo{i-1} = legenda{i-1};
%     legend(legendInfo)
%     hold on
% end

% subplot(3,3,1),plot(t,x(:,2))
% xlabel('Tempo')
% ylabel('Biomassa(gDW/L)')
% subplot(3,3,2),plot(t,x(:,3),'k')
% xlabel('Tempo')
% ylabel('Glicerol(g/L)')
% subplot(3,3,3),plot(t,x(:,4),'r')
% xlabel('Tempo')
% ylabel('Acetato(g/L)')
% subplot(3,3,4),plot(t,x(:,5),'g')
% xlabel('Tempo')
% ylabel('Lactato(g/L)')
% ylim([-1 1])
% subplot(3,3,5),plot(t,x(:,7),'m')
% xlabel('Tempo')
% ylabel('Plasmídeo(g/L)')
% subplot(3,3,6),plot(t,x(:,8),'y')
% xlabel('Tempo')
% ylabel('Asparaginase(g/L)') 
%  
% for i=1:(length(x0)-2)
%     figure(i)
%     plot(t,x(:,i))
%     xlabel('Tempo')
%     switch i
%         case 1
%             ylabel('Volume(L)')
%         case 2
%             ylabel('Biomassa(gDW/L)')
%         case 3
%             ylabel('Glicerol(g/L)')
%         case 4
%             ylabel('Acetato(g/L)')
%         case 5
%             ylabel('Lactato(g/L)')
% %             ylim([-1 1])
%         case 6
%             ylabel('Lactose(g/L)')
% %             ylim([-1 1])
%         case 7
%             ylabel('Plasmídeo(g/L)')
%         case 8
%             ylabel('Asparaginase(g/L)')    
% %         case 9
% %             ylabel('Oxigênio na fase líquida(g/L)')   
% %         case 10
% %             ylabel('Oxigênio na fase gasosa(g/L)')   
%     end
%     
% end
% 
% n_o2=length(x0)-1;
% figure(n_o2)
% plot(t,p_o2)
% xlabel('Tempo')
% ylabel('pO2')

figure(3)
subplot(2,4,1),plot(t,x(:,2))
xlabel('Tempo')
ylabel('Biomassa(gDW/L)')
subplot(2,4,2),plot(t,x(:,3))
xlabel('Tempo')
ylabel('Glicerol(g/L)')
subplot(2,4,3),plot(t,x(:,4))
xlabel('Tempo')
ylabel('Acetato(g/L)')
subplot(2,4,4),plot(t,x(:,5))
xlabel('Tempo')
ylabel('Lactato(g/L)')
subplot(2,4,5),plot(t,x(:,6))
xlabel('Tempo')
ylabel('Lactose(g/L)')
subplot(2,4,6),plot(t,x(:,7))
xlabel('Tempo')
ylabel('Plasmídeo(g/L)')
subplot(2,4,7),plot(t,Enz)
xlabel('Tempo')
ylabel('Asparaginase(UI/mL)')    
subplot(2,4,8),plot(t,p_o2)
xlabel('Tempo')
ylabel('po2')
