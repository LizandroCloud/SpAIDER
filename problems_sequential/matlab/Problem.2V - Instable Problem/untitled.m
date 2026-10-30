% Dados de exemplo
X = randn(100, 38);  % 100 amostras, 38 variáveis de estado

% 1. Centralização dos dados
X_centered = X - mean(X, 1);

% 2. Cálculo da matriz de covariância
C = (X_centered' * X_centered) / (size(X, 1) - 1);

% 3. Decomposição em valores e vetores próprios
[U, S, V] = svd(C);

% 4. Seleção dos primeiros k modos principais
k = 5;  % Número de modos principais a serem mantidos
U_k = U(:, 1:k);
S_k = S(1:k, 1:k);

% 5. Projeção dos dados no espaço reduzido
X_reduced = X_centered * U_k;

% 6. Reconstrução dos dados aproximados
X_reconstructed = X_reduced * U_k';

% 7. Visualização da energia explicada por cada modo
energy_explained = diag(S) / sum(diag(S));
disp('Energia explicada por cada modo:');
disp(energy_explained);
% Dados de exemplo
X = randn(100, 38);  % 100 amostras, 38 variáveis de estado

% 1. Centralização dos dados
X_centered = X - mean(X, 1);

% 2. Cálculo da matriz de covariância
C = (X_centered' * X_centered) / (size(X, 1) - 1);

% 3. Decomposição em valores e vetores próprios
[U, S, V] = svd(C);

% 4. Seleção dos primeiros k modos principais
k = 5;  % Número de modos principais a serem mantidos
U_k = U(:, 1:k);
S_k = S(1:k, 1:k);

% 5. Projeção dos dados no espaço reduzido
X_reduced = X_centered * U_k;

% 6. Reconstrução dos dados aproximados
X_reconstructed = X_reduced * U_k';

% 7. Visualização da energia explicada por cada modo
energy_explained = diag(S) / sum(diag(S));
disp('Energia explicada por cada modo:');
disp(energy_explained);
