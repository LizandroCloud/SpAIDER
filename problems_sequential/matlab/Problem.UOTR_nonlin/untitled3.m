% Dados de exemplo
z = iddata(yx, ux, 0.2);
ze = z(1:300, :, 1);
zv = z(301:541, :, 1);
n_inputs = size(ux, 2);  % Número de entradas
n_outputs = size(yx, 2); % Número de saídas

% Criar o modelo NLARX
% m4 = idnlarx([2, 2, n_outputs]);
% Criação do modelo
% m4 = idnlarx([2 2 3], idSigmoidNetwork(14));
m4 = nlarx(ze, nlarx1);  % Treinamento

% Inicialize um vetor de predições
n_samples = size(zv, 1);
n_outputs = size(zv.y, 2);  % Número de saídas
predictions = zeros(n_samples, n_outputs);

% Realize a predição
for i = 1:n_samples
    input_vector = zv.u(i, :);  % Dados de entrada correspondentes
    predictions(i, :) = predict(m4, input_vector);
end

% Visualize os resultados
figure;
hold on;
for j = 1:n_outputs
    subplot(n_outputs, 1, j);
    plot(zv.y(:, j), 'b--', 'DisplayName', 'Valores Reais');
    hold on;
    plot(predictions(:, j), 'r', 'DisplayName', 'Predições');
    xlabel('Amostra');
    ylabel(['Saída ' num2str(j)]);
    legend;
end
sgtitle('Predições com Modelo NLARX');
