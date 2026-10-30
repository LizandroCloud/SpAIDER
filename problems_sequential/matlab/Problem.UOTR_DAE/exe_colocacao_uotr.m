

global PP
PP=blackoil;
global  modelo_x0 
load('uotr_ctrl_ss_v1.mat','x0','u0','par');


% idu=[1,5,7,9,15,16,19,20,25];
idu=[1];
t0 = 0; % Initial time for integration
tf = 3600*2; % Final time for integration
x0 = modelo_x0; %state Initial Condition


% Definir intervalo de tempo e condição inicial
tspan = [0, tf];
y0 = x0;  % Condição inicial para 4 EDOs
h = 0.01;

% Chamar método de colocação de Gauss-Legendre
tic
[T, Y] = gauss_legendre_collocation(@systemODE, tspan, y0, h);
toc

% Chamar ode45 para resolver o sistema de EDOs
tic
[t, y] = ode15s(@systemODE, tspan, y0);
toc
% Plotar resultados
figure;
hold on;
for i = 1:size(Y, 2)
    plot(T, Y(:, i), 'DisplayName', ['y' num2str(i)]);
end
xlabel('Tempo');
ylabel('Solução');
legend;
title('Solução do Sistema de 4 EDOs usando o Método de Colocação de Gauss-Legendre');

% Plotar resultados
figure;
hold on;
plot(t, y(:, 1), 'DisplayName', 'y1');
plot(t, y(:, 2), 'DisplayName', 'y2');
plot(t, y(:, 3), 'DisplayName', 'y3');
plot(t, y(:, 4), 'DisplayName', 'y4');
xlabel('Tempo');
ylabel('Solução');
legend;
title('Solução do Sistema de 4 EDOs usando ode45');


function dydt = systemODE(t, y)
    % Exemplo de sistema de 4 EDOs
    global modelo_u0 modelo_par

    ws_reg=ver_reg_est1(ws,yi,t_out);

    %idx=[2,3,5,6,7,9,10,13,15,17,19,20,22,24,25,28];
    %idy=[1,4,8,12,14,16,18,26]; %idy=[1,4,8,11,12,14,16,18,21,23,26,27];

    idu=[1]; %idu=[1,5,7,9,14,15,16,19,20,21,23,25,26];
    %idw=[2,3,4,6,8,10,11,12,13,17,18,22,24,27];

    u=modelo_u0; u(idu)=ws_reg*1e3;
    dx=uotr_ctrl(t_out,x,u,1,[],modelo_par);
    dx=dx*tal;
end


function [T, Y] = gauss_legendre_collocation(odefun, tspan, y0, h)
    % tspan: intervalo de tempo [t0, tf]
    % y0: condição inicial (vetor de tamanho N)
    % h: tamanho do passo

    t0 = tspan(1);
    tf = tspan(2);
    N = (tf - t0) / h;  % Número de passos
    
    T = t0:h:tf;
    M = length(y0);  % Número de EDOs
    Y = zeros(length(T), M);  % Matriz de soluções
    Y(1, :) = y0;

    % Pontos e pesos de Gauss-Legendre para 2 pontos
    c1 = 1 / 2;
    c2 = 1 / 2;
    w1 = 1;
    w2 = 1;

    for n = 1:N
        tn = T(n);
        yn = Y(n, :)';  % Transpor para tornar yn um vetor coluna

        % Pontos de colocação de Gauss-Legendre
        k1 = odefun(tn + c1 * h, yn);
        k2 = odefun(tn + c2 * h, yn + h * k1);
        
        % Atualizar solução
        Y(n+1, :) = (yn + h * (w1 * k1 + w2 * k2))';  % Transpor para tornar um vetor linha
    end
end
