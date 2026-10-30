% Definir o sistema de EDO
% Condição inicial
y0 = 0;

% Intervalo de tempo
tspan = [0 10];

% Pontos específicos para uma malha fixa
t_fixed = linspace(0, 10, 100);  % Exemplo de 100 pontos

% Resolver a EDO com uma malha fixa
tic;
[t_fixed_sol, y_fixed] = ode15s(@odefun, t_fixed, y0);
toc;

% Resolver a EDO adaptativamente
tic;
[t_adaptive, y_adaptive] = ode15s(@odefun, tspan, y0);
toc;

% Comparar os resultados
disp('Resultados com malha fixa:');
disp(['Tempo total: ', num2str(toc)]);

disp('Resultados com método adaptativo:');
disp(['Tempo total: ', num2str(toc)]);

% Interpolação para comparação
y_interpolated = interp1(t_adaptive, y_adaptive, t_fixed, 'linear');

function dydt = odefun(t, y)
    dydt = -2*y + 1;
end
