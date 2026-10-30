% Suponha que 'nlarx_model' é seu modelo NLARX já criado
% E 'input_data' e 'output_data' são seus dados de entrada e saída

% Defina os dados de entrada para a predição
% input_data deve ser um vetor ou matriz com os dados exógenos
% output_data deve conter os dados anteriores, se necessário
input_data = ux;
output_data = yx;
nlarx_model = nlarx1;

% Supondo que nlarx_model é seu modelo já treinado

% Dados de entrada (X) - matriz com dimensão [n_samples, n_inputs]
% Dados de saída (Y) - vetor ou matriz com dimensão [n_samples, n_outputs]

% Verifique as dimensões
[n_samples, n_inputs] = size(input_data);
[n_samples_out, n_outputs] = size(output_data);

n = length(input_data);  % Número de amostras

% Previsão para o horizonte de previsão desejado
h = 10;  % Horizonte de previsão
predictions = zeros(h, 38);  % Inicialize o vetor de predições

% Realize a predição
for i = 1:h
    % Para a primeira predição, use o último dado conhecido de input_data
    if i == 1
        input_vector = input_data(end);  % Último valor do input
    else
        % Para previsões subsequentes, você pode usar o último valor de saída previsto
        input_vector = input_data(end);  % Último dado de entrada (continua o mesmo)
    end
    
    % Predição
    predictions(i, :) = predict(nlarx_model, input_vector);
    
    % Atualize os dados de entrada (opcional)
    % Se o modelo exige a adição de previsões anteriores como entradas, você deve adaptar isso
    input_data = [input_data; input_vector];  % Para continuar usando se necessário
end

% Visualize os resultados
figure;
hold on;
for j = 1:38
    subplot(38, 1, j);
    plot(1:h, predictions(:, j), 'r', 'DisplayName', 'Predição');
    hold on;
    plot(1:h, output_data(end-h+1:end, j), 'b--', 'DisplayName', 'Valores Reais');
    xlabel('Passos de Tempo');
    ylabel(['Saída ' num2str(j)]);
    legend;
end
sgtitle('Predições com Modelo NLARX');