function [ pF, vF, p ] = filter_pox(tss, ww, point)
p=point;
% Vetores de entrada
% valores = [1, 2, 2, 2, 2, 3, 4, 4, 4, 5];
% posicoes = 1:length(valores);  % ou outro vetor de posições correspondente
valores = round(ww,2);
posicoes = round(tss,2);
% Inicializar vetor lógico para manter todos os índices
indicesParaManter = true(size(valores));

i = 1;
while i <= length(valores)
    j = i;
    
    % Procurar sequência de valores consecutivos iguais
    while j < length(valores) && valores(j) == valores(j + 1)
        j = j + 1;
    end
    
    comprimento = j - i + 1;
    
    if comprimento >= 3
        % Remove apenas os elementos intermediários (mantém os extremos)
        indicesParaManter((i+1):(j-1)) = false;
        
    end
    
    i = j + 1;
end

% Aplicar a filtragem aos dois vetores
vF = valores(indicesParaManter);
pF = posicoes(indicesParaManter);
p = point(indicesParaManter);

% % Mostrar resultados
% disp('Valores filtrados:');
% disp(valoresFiltrados);
% disp('Posições filtradas:');
% disp(posicoesFiltradas);
