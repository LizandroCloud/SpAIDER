% Gerar dados fictícios
n_samples = 100; % Número de amostras
n_outputs = 1;  % Número de saídas
n_inputs = 1;   % Número de entradas

% Criar dados de entrada (ux) e saída (yx)
ux = rand(n_samples, n_inputs); % Dados de entrada aleatórios
yx = rand(n_samples, n_outputs); % Dados de saída aleatórios

% Criar o objeto iddata
z = iddata(yx, ux, 1); % 1 segundo de intervalo de amostragem

% Dividir os dados em conjunto de treinamento e validação
ze = z(1:70, :, 1);  % 7 amostras para treinamento
zv = z(71:100, :, 1); % 3 amostras para validação

% Verifique as dimensões
disp('Dimensões de ze.u:');
disp(size(ze.u));  % Dimensões das entradas
disp('Dimensões de ze.y:');
disp(size(ze.y));  % Dimensões das saídas

% Definir o modelo NLARX
% p = ordem de saída, q = ordem de entrada, r = número de saídas
p = 2; % Ordem da saída
q = 1; % Ordem da entrada
r = n_outputs; % Número de saídas
m4 = idnlarx([p, q, r], idSigmoidNetwork(10));  % Criação do modelo

% Treinar o modelo
try
    m4 = nlarx(ze, m4);  % Treinamento
catch ME
    disp('Erro ao treinar o modelo:');
    disp(ME.message);
end

% Realizar previsões
n_samples_val = size(zv, 1);
predictions = zeros(n_samples_val, n_outputs);

for i = 1:n_samples_val
    input_vector = zv.u(i, :); % Pegar os dados de entrada
    predictions(i, :) = predict(m4, input_vector'); % Prever
end

% Visualizar os resultados
figure;
hold on;
for j = 1:n_outputs
    subplot(n_outputs, 1, j);
    plot(zv.y(:, j), 'b--', 'DisplayName', 'Valores Reais'); % Valores reais
    hold on;
    plot(predictions(:, j), 'r', 'DisplayName', 'Predições'); % Predições
    xlabel('Amostra');
    ylabel(['Saída ' num2str(j)]);
    legend;
end
sgtitle('Predições com Modelo NLARX');
