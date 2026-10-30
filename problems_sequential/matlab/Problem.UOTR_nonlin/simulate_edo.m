% Função para simular a EDO
function [t, x] = simulate_edo(u, tspan)
    x0 = [0; 0]; 
    ode_func = @(t, x) [-x(1) + u; -x(2) + x(1)];
    [t, x] = ode45(@(t, x) ode_func(t, x), tspan, x0);
end

% Principal
tspan = [0 10]; 
u = @(t) sin(t); % Exemplo de entrada
[t, x] = simulate_edo(@(t) u(t), tspan);

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
