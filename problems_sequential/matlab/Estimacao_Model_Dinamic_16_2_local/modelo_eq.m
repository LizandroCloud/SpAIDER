function [ dx ] = modelo_eq( t_out, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg)
[u] = ver_reg_est1(ws,yi,t_out);
[uref] = ver_reg_est1(ws,yi,t_out);
%global epsilon V_g V_L F_ar c_m c_g0 t_i t_f

epsilon=1e-6;
T=37;
R=0.082; %atm L/mol K
P=1; %atm
x_o2=0.21;
V_g=1;%L
V_L=1.5;%L
F_ar=90;%L/h
t_i = 8;


k_ref=0.8;
kd_ref=0.01;
alfa_ref=50;
km_glyc_ref=0.05;
km_lact_ref=0.05;
vmax_glyc_ref=10;
vmax_lact_ref=1;
v_max_o2_ref=20000;
km_o2_ref=0.0024;
kla_ref=41.3/5;
v_lactato_ref=5;

k=k_ref*u(1);
kd=kd_ref*u(2);
alfa=alfa_ref*u(3);
km_glyc=km_glyc_ref*u(4);
km_lact=km_lact_ref*u(5);
vmax_glyc=vmax_glyc_ref*u(6);
vmax_lact=vmax_lact_ref*u(7);
v_max_o2=v_max_o2_ref*u(8);
km_o2=km_o2_ref*u(9);
kla=kla_ref*u(10);
v_lactato=v_lactato_ref*u(11);

load('iJO1366_PM_asp_plasm_novo_2.mat');
modelo_1=iJO1366_PM_asp_plasm_novo_2;

modelo_1 = changeRxnBounds(modelo_1,'EX_glc(e)',0,'b');
modelo_1 = changeRxnBounds(modelo_1,'EX_glyc(e)',(-7),'l'); 
modelo_1 = changeRxnBounds(modelo_1,'EX_glyc(e)',0,'u'); 
modelo_1 = changeRxnBounds(modelo_1,'NADTRHD',4,'u'); %NAD transhydrogenase
modelo = changeRxnBounds(modelo_1,'G6PDH2r',3.5,'u');


modelo=constraints(x,modelo,km_glyc,km_lact,vmax_glyc,vmax_lact,v_max_o2,km_o2,v_lactato);

ID_o2 = findMetIDs(modelo,'o2[e]'); 

% Conc. de o2 na vazão de alimentação de ar (c_g0)
c_g0=P*modelo.PM(ID_o2)*x_o2/(R*(T+273.15));

% Conc. de saturação de o2 na fase líquida

c_m=3.62e-06*T^2-3.312e-04*T+0.01466217;



if t_out>t_i && t_out<=t_f
   f_feed=10; % L/h
elseif t_out<=t_i
    f_feed=0;
elseif t_out>t_f
    f_feed=0;
end

[v_class,obj,erro]=fba2(modelo,[],[],alfa,1);

[v,obj2,erro]=fba2_multi(modelo,k,obj,alfa,epsilon,2);

mi=v(8);

% [ID_glyc, v_glyc, ID_ac, v_ac, ID_lac, v_lac,ID_lactose,v_lactose]=fluxos(modelo,v);
[ID,f]=fluxos(modelo,v); % onde f é o fluxo

G_feed=0;%conc. de alimentação de glicerol
% f_feed=0;

% Equações diferenciais
    
dx(1)=f_feed; %dV/dt - [L/h]

dx(2)= (mi-kd)*x(2) - x(2)*f_feed/x(1); %Biomassa dX/dt - [gDW/L h]

dx(3)= f_feed*(G_feed-x(3))/x(1) + f(1)*modelo.PM(ID(1))*x(2)/1000; % Substrato: glicerol - [g/L h]

dx(4)= f(2)*modelo.PM(ID(2))*x(2)/1000 - x(4)*f_feed/x(1); % Acetato - [g/L h]

dx(5)= f(3)*modelo.PM(ID(3))*x(2)/1000 - x(5)*f_feed/x(1); % Lactato - [g/L h]

% %variável binária relacionada a lactose
%  if f(4)~=0 %v_lactose
%      b=1;
%  else
%      b=0;
%  end

% b=1;
% L_feed=0;
L_feed=75;%g/L

%Lactose-D
dx(6)= f_feed*(L_feed-x(6))/x(1) + f(4)*modelo.PM(ID(4))*x(2)/1000; 

%Plasmídeo
dx(7)=f(6)*modelo.PM(ID(6))*x(2)/1000 - x(7)*f_feed/x(1);

% Asparaginase
% dx(8)=b*(f(5)*modelo.PM(ID(5))*x(2)/1000 - x(8)*f_feed/x(1));

dx(8)=f(5)*modelo.PM(ID(5))*x(2)/1000 - x(8)*f_feed/x(1); %[g/L h]

%Considerando a seguinte conversão: 1mg da enzima - 300 UI de atividade
%enzimática
% Fator de conversão: 300 (chute inicial)
% Unidade final: UI/mL

%Constante de equlíbrio
K=c_m/c_g0;
% Conc. de o2 na fase líquida em equilíbrio com o gás (c_le)
c_le=K*x(10);

% O2 na fase líquida:
dx(9)=kla*(c_le-x(9))-f(7)*modelo.PM(ID(7))*x(2)/1000;

% O2 na fase gasosa:
dx(10)=F_ar*(c_g0-x(10))/V_g-kla*(c_le-x(9))*V_L/V_g;

dx = dx*tal;

dx=transp(dx);

end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
