% Exemplo de uso
A = [1, 2; 3, 4]; % Matriz A (2x2)
d = 2; % Dimensão adicional
B = nonlinear_mapping(A, d);
disp('Matriz B resultante:');
disp(B);

function B = nonlinear_mapping(A, d)
    % A é uma matriz n x n
    % d é a dimensão adicional
    
    [n, ~] = size(A);
    
    % Inicializa a matriz B (n+d) x (n+d)
    B = zeros(n + d, n + d);
    
    % Mapeamento não linear (exemplo: polinômio quadrático)
    for i = 1:n
        for j = 1:n
            % Preservar a posição original
            B(i, j) = A(i, j);
            
            % Adicionar novas características
            B(i, n + 1) = A(i, j)^2;  % Exemplo de mapeamento quadrático
            B(n + 1, j) = A(i, j)^2;  % Preservando a simetria
            % Adicione mais mapeamentos conforme necessário
        end
    end
    
    % Adicionar mais linhas e colunas conforme necessário para d
    for k = 1:d
        B(n + k, n + k) = 1; % Adicionando valores na diagonal (exemplo)
    end
end


