% Definir intervalo de tempo e condição inicial
tspan = [0, 5];
y0 = [1; 1; 1; 1];  % Condição inicial para 4 EDOs
h = 0.1;

% Chamar método de colocação de Gauss-Legendre
[T, Y] = gauss_legendre_collocation(@systemODE, tspan, y0, h);

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

function dydt = systemODE(t, y)
    % Exemplo de sistema de 4 EDOs
    dydt = zeros(4, 1);
    dydt(1) = -2 * y(1) + t^2;
    dydt(2) = -3 * y(2) + sin(t);
    dydt(3) = -y(3) + cos(t);
    dydt(4) = y(4)^2 - t;
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
        yn = Y(n, :)';

        % Pontos de colocação de Gauss-Legendre
        k1 = odefun(tn + c1 * h, yn);
        k2 = odefun(tn + c2 * h, yn + h * k1);
        
        % Atualizar solução
        Y(n+1, :) = yn + h * (w1 * k1 + w2 * k2)';
    end
end


