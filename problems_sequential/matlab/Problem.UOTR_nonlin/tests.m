global PP
PP=blackoil;
global  modelo_x0 
load('uotr_ctrl_ss_v1.mat','x0','u0','par');


% idu=[1,5,7,9,15,16,19,20,25];
idu=[1];
t0 = 0; % Initial time for integration
tf = 7200; % Final time for integration
x0 = modelo_x0; %state Initial Condition


% Definir intervalo de tempo e condição inicial
tspan = [0, tf];
y0 = x0;  % Condição inicial para 4 EDOs

tic
u0=18; tx=[]; yx=[]; tend=0; ux0=[];
for k=1:10
    u=u0+5*(-1)^(k);
    [t, y] = ode15s(@systemODE, tspan/10, y0, [], u);
    tx=[tx;t(1:end)+tend];
    yx=[yx;y(1:end,:)];
    ux=[ux0; ones(length(t),1)*u];
    ux0=ux;
    y0=y(end,:);
    tend = tx(end);
end
toc
figure;
hold on;
for i = 1:size(yx, 2)
    plot(tx, yx(:, i), 'DisplayName', ['y' num2str(i)]);
end

xlabel('Tempo');
ylabel('Solução');
legend;
title('Solução do Sistema de 4 EDOs usando ode45');




function dydt = systemODE(t_out, x, uz)
    % Exemplo de sistema de 4 EDOs
    global modelo_u0 modelo_par

    % ws_reg=ver_reg_est1(ws,yi,t_out);

    %idx=[2,3,5,6,7,9,10,13,15,17,19,20,22,24,25,28];
    %idy=[1,4,8,12,14,16,18,26]; %idy=[1,4,8,11,12,14,16,18,21,23,26,27];

    idu=[1]; %idu=[1,5,7,9,14,15,16,19,20,21,23,25,26];
    %idw=[2,3,4,6,8,10,11,12,13,17,18,22,24,27];

    u=modelo_u0; u(idu)=uz*1e3;
    dx=uotr_ctrl(t_out,x,u,1,[],modelo_par);
    dydt=dx*1;
end