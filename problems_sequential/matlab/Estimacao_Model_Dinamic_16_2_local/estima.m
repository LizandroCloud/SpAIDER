function estima

% Estimacao de parametros 

global x0 M1 M2 tempo modelo t_i t_f V_g V_L F_ar c_m c_g0 epsilon contador

contador=0;

t_i = 8;
t_f = 8.01;

load('iJO1366_PM_asp_plasm_novo_2.mat');
modelo_1=iJO1366_PM_asp_plasm_novo_2;

modelo_1 = changeRxnBounds(modelo_1,'EX_glc(e)',0,'b');
modelo_1 = changeRxnBounds(modelo_1,'EX_glyc(e)',(-7),'l'); 
modelo_1 = changeRxnBounds(modelo_1,'EX_glyc(e)',0,'u'); 
modelo_1 = changeRxnBounds(modelo_1,'NADTRHD',4,'u'); %NAD transhydrogenase
modelo = changeRxnBounds(modelo_1,'G6PDH2r',3.5,'u');

epsilon=1e-6;
% kla=41.3;
T=37;
R=0.082; %atm L/mol K
P=1; %atm
x_o2=0.21;
V_g=1;%L
V_L=1.5;%L
F_ar=90;%L/h

ID_o2 = findMetIDs(modelo,'o2[e]'); 

% Conc. de o2 na vazão de alimentação de ar (c_g0)
c_g0=P*modelo.PM(ID_o2)*x_o2/(R*(T+273.15));

% Conc. de saturação de o2 na fase líquida
c_m=3.62e-06*T^2-3.312e-04*T+0.01466217;


V0=1.5;
X0=5.3640E-02;
G0=3.1337E+01;
Ac0=0;
Lac0=0;
L0=0;
Asp0=0;
Plasm0=0;
C_L0=0.0072;
C_G0=c_g0;

x0=[V0 X0 G0 Ac0 Lac0 L0 Plasm0 Asp0 C_L0 C_G0];


% k1=1;
% kd=1;
% alfa=1;
% km_glyc=1;
% km_lact=1;
% vmax_glyc=1;
% vmax_lact=1;
% v_max_o2=1;
% km_o2=1;
% kla=1;
% v_lactato=1;
% 
% par0=[k1 kd alfa km_glyc km_lact vmax_glyc vmax_lact v_max_o2 km_o2 kla v_lactato];

load par_inicial.mat

par0=par_inicial;

np=length(par0);

biomassa=[
        0.0000E+00     5.3640E-02
        2.0000E+00     9.9617E-02
        4.0000E+00     2.6820E-01
        6.0000E+00     7.9693E-01
        7.0000E+00     1.1801E+00
        8.0000E+00     1.3027E+00
        1.1000E+01     1.6935E+00
        1.4000E+01     1.8467E+00
        1.7000E+01     1.7471E+00
        2.0000E+01     1.7241E+00
        2.3000E+01     1.6245E+00];

glicerol=[
    0.0000E+00     3.1337E+01
    2.0000E+00     3.1744E+01
    4.0000E+00     3.0659E+01
    6.0000E+00     2.9438E+01
    8.0000E+00     2.4554E+01
    1.1000E+01     2.0756E+01
    1.4000E+01     1.6143E+01
    1.7000E+01     1.3430E+01
    2.0000E+01     1.0988E+01
    2.3000E+01     9.3605E+00];

lactose=[
    0.0000E+00     0.0000E+00
    2.0000E+00     0.0000E+00
    4.0000E+00     0.0000E+00
    6.0000E+00     0.0000E+00
    8.0000E+00     0.0000E+00
    8.0100E+00     4.4574E+00
    1.1000E+01     4.3217E+00
    1.4000E+01     3.4302E+00
    1.7000E+01     2.2674E+00
    2.0000E+01     8.5271E-01
    2.3000E+01     9.6899E-02
  ];

acetato=[                                            
    0.0000E+00     0.0000E+00
    2.0000E+00     0.0000E+00
    4.0000E+00     0.0000E+00
    6.0000E+00     0.0000E+00
    8.0000E+00     4.9080E-01
    1.1000E+01     9.0491E-01
    1.4000E+01     9.2025E-01
    1.7000E+01     1.1656E+00
    2.0000E+01     1.4417E+00
    2.3000E+01     1.6718E+00];

l_asp=[
    0.0000E+00     0.0000E+00
    2.0000E+00     0.0000E+00
    4.0000E+00     0.0000E+00
    6.0000E+00     0.0000E+00
    8.0000E+00     0.0000E+00
    1.1000E+01     0.0000E+00
    1.4000E+01     2.6336E+00
    1.7000E+01     6.2977E+00
    2.0000E+01     8.6641E+00
    2.3000E+01     7.1374E+00];

o2=[0.0000E+00     9.7719E+01
    2.0000E+00     9.8099E+01
    4.0000E+00     6.5779E+01
    6.0000E+00     1.3688E+01
    7.0000E+00     1.4829E+01
    8.0000E+00     2.2814E+00
    1.1000E+01     2.6616E+00
    1.4000E+01     3.0418E+00
    1.7000E+01     3.0418E+00
    2.0000E+01     3.8023E+00
    2.3000E+01     2.6616E+00];

lactato=[0.0000E+00    0.0000E+00
        2.0000E+00     0.0000E+00
        4.0000E+00     0.0000E+00
        6.0000E+00     0.0000E+00
        8.0000E+00     0.0000E+00
        1.1000E+01     2.9141E+00
        1.4000E+01     6.8179E+00
        1.7000E+01     1.0447E+01
        2.0000E+01     1.2976E+01
        2.3000E+01     1.4296E+01];

tempo=0:1:23;
M1={biomassa(:,1) glicerol(:,1) lactose(:,1) acetato(:,1) o2(:,1) l_asp(:,1) lactato(:,1)};
M2={biomassa(:,2) glicerol(:,2) lactose(:,2) acetato(:,2) o2(:,2) l_asp(:,2) lactato(:,2)};

% options = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
%  'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 1e-5,...
%  'TolFun', 1e-5, 'TolCon', 1e-5, 'MaxFunEval', 15000000, 'MaxIter', 1e6);
% lb=ones(1,np)*1e-4;
% ub=ones(1,np)*10;
% ub(1)=1;
% paropt = fmincon('minimiza', par0, [], [], [], [],lb, ub, [] , options);  

% paropt=buscarndb('minimiza',par0,0,1,100,ones(1,np)*1e-4,ones(1,np)*20,[],[],[],2*ones(1,np));
%paropt=swarm('minimiza',par0,0,ones(1,np),ones(1,np)*1e-4,ones(1,np)*20,2*ones(1,np),100,20,[],1,1);
% options = optimoptions('particleswarm','SwarmSize',50,'HybridFcn',@fmincon,'Display','iter','ObjectiveLimit',100);

% [paropt,f_obj,exitflag,output] = particleswarm(@(par) minimiza(par),np,lb,ub,options);

load('paropt.mat')

% save('paropt.mat','paropt')

op=odeset('MaxStep',0.01);
op1=odeset('AbsTol',1e-5);
[t,x]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,paropt),[0:1:t_i],x0);
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,paropt),[t_i:0.01:t_f],x(end,:),op);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,paropt),[t_f:0.00001:11],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,paropt),[11:0.00001:15],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,paropt),[15:0.00001:18],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,paropt),[18:0.00001:20],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,paropt),[20:0.00001:21],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,paropt),[21:0.00001:22],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,paropt),[22:0.00001:23],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 
p_o2=100*x(:,9)/c_m;
Enz=pseudo_steady_state_v(x);

figure(1)
plot(biomassa(:,1),biomassa(:,2),'ob',t,x(:,2),'-r')
xlabel('Tempo (h)')
ylabel('Biomassa')
figure(2)
plot(glicerol(:,1),glicerol(:,2),'ob',t,x(:,3),'-r')
xlabel('Tempo (h)')
ylabel('Glicerol')
figure(3)
plot(lactose(:,1),lactose(:,2),'ob',t,x(:,6),'-r')
xlabel('Tempo (h)')
ylabel('Lactose')
figure(4)
plot(acetato(:,1),acetato(:,2),'ob',t,x(:,4),'-r')
xlabel('Tempo (h)')
ylabel('Acetato')
figure(5)
plot(l_asp(:,1),l_asp(:,2),'ob',t,Enz,'-r')
xlabel('Tempo (h)')
ylabel('L-asparaginase')
figure(6)
plot(o2(:,1),o2(:,2),'ob',t,p_o2,'-r')
xlabel('Tempo (h)')
ylabel('pO2')
figure(7)
plot(lactato(:,1),lactato(:,2),'ob',t,x(:,5),'-r')
xlabel('Tempo (h)')
ylabel('Lactato')