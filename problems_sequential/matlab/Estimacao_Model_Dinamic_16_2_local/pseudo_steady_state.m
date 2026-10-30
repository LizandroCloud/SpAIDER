% Modelo dinâmico - Pseudo-estacionário
function dx=pseudo_steady_state(t,x,modelo,par)

global epsilon V_g V_L F_ar c_m c_g0 t_i t_f

k_ref=0.8;
kd_ref=0.1;
alfa_ref=50;
km_glyc_ref=0.05;
km_lact_ref=0.05;
vmax_glyc_ref=10;
vmax_lact_ref=1;
v_max_o2_ref=20000;
km_o2_ref=0.0024;
kla_ref=41.3/5;
v_lactato_ref=5;

k=k_ref*par(1);
kd=kd_ref*par(2);
alfa=alfa_ref*par(3);
km_glyc=km_glyc_ref*par(4);
km_lact=km_lact_ref*par(5);
vmax_glyc=vmax_glyc_ref*par(6);
vmax_lact=vmax_lact_ref*par(7);
v_max_o2=v_max_o2_ref*par(8);
km_o2=km_o2_ref*par(9);
kla=kla_ref*par(10);
v_lactato=v_lactato_ref*par(11);

modelo=constraints(x,modelo,km_glyc,km_lact,vmax_glyc,vmax_lact,v_max_o2,km_o2,v_lactato);

if t>t_i && t<=t_f
   f_feed=10; % L/h
elseif t<=t_i
    f_feed=0;
elseif t>t_f
    f_feed=0;
end

t

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

dx=dx';
