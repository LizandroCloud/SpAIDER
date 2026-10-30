% Função para simular a EDO

% Principal
tspan = [0 10]; 
% u = @(t) sin(t); % Exemplo de entrada
[t, x] = ode45(@(t, x) ode_func(t, x), tspan, x0);

% Saída baseada nos estados
y = x(:, 1) + 2 * x(:, 2); 

% Criar conjunto de dados
data = iddata(y, u(t), 0.1); % 0.1 é o intervalo de amostragem

% Criar o modelo NLARX
na = 2; % Ordem da saída
nb = 2; % Ordem da entrada
nk = 0; % Atraso
model = nlarx(data, 'nlarx', [na nb nk]);

% Validação do modelo
compare(data, model);

function [dx] = ode_func(t,x)
    u= sin(t);
    dx(1) = -x(1) + u;
    dx(2) = -x(2) + x(1);
end

function [t, x] = simulate_edo(tspan)
    x0 = [0; 0]; 
    % ode_func = @(t, x) [-x(1) + u; -x(2) + x(1)];
    [t, x] = ode45(@(t, x) ode_func(t, x), tspan, x0);
end