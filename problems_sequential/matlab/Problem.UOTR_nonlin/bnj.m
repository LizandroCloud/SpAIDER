% Exemplo de uso
A = [1, 2; 3, 4]; % Matriz A (2x2)
d = 2; % Dimensão adicional
B = interpolate_symmetric_mapping(A, d);
disp('Matriz B resultante:');
disp(B);

function B = interpolate_symmetric_mapping(A, d)
    % A é uma matriz simétrica n x n
    [n, ~] = size(A);
    
    % Dimensões da nova matriz
    new_size = n + d;
    B = zeros(new_size, new_size);
    
    % Criar um grid de interpolação
    [X, Y] = meshgrid(1:n, 1:n); % Coordenadas originais
    [Xq, Yq] = meshgrid(linspace(1, n, new_size), linspace(1, n, new_size)); % Nova grade
    
    % Interpolar valores
    B = interp2(X, Y, A, Xq, Yq, 'linear'); % Interpolação linear
    
    % Garantir a simetria
    B = (B + B') / 2; % Força a matriz a ser simétrica
end

