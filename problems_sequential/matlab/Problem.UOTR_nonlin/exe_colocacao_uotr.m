

global PP
PP=blackoil;
global  modelo_x0 
load('uotr_ctrl_ss_v1.mat','x0','u0','par');


% idu=[1,5,7,9,15,16,19,20,25];
idu=[1];
t0 = 0; % Initial time for integration
tf = 3600; % Final time for integration
x0 = modelo_x0; %state Initial Condition


% Definir intervalo de tempo e condição inicial
tspan = [0, tf];
tspan = linspace(0,tf,10);
y0 = x0;  % Condição inicial para 4 EDOs


% Chamar método de colocação de Gauss-Legendre
Z=[]; nd=10;
h = tspan(end)/nd;
t=[0];
tic
tol=1e-3;
uz=24e3;
for i=1:nd
    
% [T, Y] = gauss_legendre_collocation(@systemODE, tspan/nd, y0, h);
[T, Y] = quadrature_collocation(@systemODE, tspan(end)/nd, y0, h, uz);
% y0=Y(end,:);
Z = [Z; Y];
S = T+t(end);
t = [t S];
end
toc
Y=Z;
T=t(2:end);
% Chamar ode45 para resolver o sistema de EDOs
tic
u=24e3;
[t, y] = ode15s(@systemODE, tspan, x0, [], u);
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
% plot(t, y(:, 1), 'DisplayName', 'y1');
% plot(t, y(:, 4), 'DisplayName', 'y2');
% plot(t, y(:, 6), 'DisplayName', 'y3');
% plot(t, y(:, 8), 'DisplayName', 'y4');
for i = 1:size(y, 2)
    plot(t, y(:, i), 'DisplayName', ['y' num2str(i)]);
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

    u=modelo_u0; u(idu)=uz;
    dx=uotr_ctrl(t_out,x,u,1,[],modelo_par);
    dydt=dx*1;
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

function [T, Y] = adaptive_rk(odefun, tspan, y0, h0, tol)
    % tspan: intervalo de tempo [t0, tf]
    % y0: condição inicial (vetor de tamanho N)
    % h0: tamanho do passo inicial
    % tol: tolerância para o controle de erro

    t0 = tspan(1);
    tf = tspan(2);
    T = t0;
    Y = y0;
    
    h = h0;
    tn = t0;
    yn = Y(:, end);

    while tn < tf
        if tn + h > tf
            h = tf - tn; % Ajustar o passo para não ultrapassar tf
        end
        
        % Estimativa usando RK4
        k1 = odefun(tn, yn);
        k2 = odefun(tn + h/2, yn + h/2 * k1);
        k3 = odefun(tn + h/2, yn + h/2 * k2);
        k4 = odefun(tn + h, yn + h * k3);
        
        y_next = yn + h/6 * (k1 + 2*k2 + 2*k3 + k4);
        
        % Estimativa com passo menor
        h_half = h / 2;
        k1_half = odefun(tn, yn);
        k2_half = odefun(tn + h_half, yn + h_half * k1_half);
        k3_half = odefun(tn + h_half, yn + h_half * k2_half);
        k4_half = odefun(tn + h_half, yn + h_half * k3_half);
        y_next_half = yn + h_half / 6 * (k1_half + 2*k2_half + 2*k3_half + k4_half);
        
        k1_half2 = odefun(tn + h_half, y_next_half);
        k2_half2 = odefun(tn + h, y_next_half + h_half * k1_half2);
        k3_half2 = odefun(tn + h, y_next_half + h_half * k2_half2);
        
        y_next_full = y_next_half + h_half / 6 * (k1_half2 + 2*k2_half2 + 2*k3_half2);
        
        % Estimar erro
        error_estimate = norm(y_next - y_next_full);
        
        if error_estimate < tol
            tn = tn + h; 
            yn = y_next;
            Y = [Y, yn]; % Armazenar a nova solução
            T = [T; tn]; % Armazenar o novo tempo
            % Aumentar o tamanho do passo
            h = h * min(2, (tol / error_estimate)^(1/4)); 
        else
            % Diminuir o tamanho do passo
            h = h * max(0.1, 0.9 * (tol / error_estimate)^(1/4));
        end
    end
    
    Y = Y'; % Transpor para formato desejado
end

function [T, Y] = quadrature_collocation(odefun, tspan, y0, h, uz)
    % tspan: intervalo de tempo [t0, tf]
    % y0: condição inicial (vetor de tamanho N)
    % h: tamanho do passo

    t0 = 0;
    tf = tspan(1);
    N = ceil((tf - t0) / h);  % Número de passos
    T = linspace(t0, tf, N + 1);  % Vetor de tempos
    M = length(y0);  % Número de EDOs
    Y = zeros(length(T), M);  % Matriz de soluções
    Y(1, :) = y0;  % Condição inicial

    % Pontos e pesos de Gauss-Legendre para 2 pontos
    c = [0.5 - 0.5/sqrt(3); 0.5 + 0.5/sqrt(3)];  % Pontos
    w = [1; 1];  % Pesos

    for n = 1:N
        tn = T(n);
        yn = Y(n, :)';  % Transpor para vetor coluna

        % Estimativas de quadratura
        k = zeros(M, 2);
        for i = 1:2
            k(:, i) = odefun(tn + c(i) * h, yn + h * sum(w(1:i)' .* k(:, 1:i), 2),uz); 
        end
        
        % Atualizar solução usando quadratura
        Y(n + 1, :) = (yn + h * sum(w' .* k, 2))';  % Transpor para vetor linha
    end
end


