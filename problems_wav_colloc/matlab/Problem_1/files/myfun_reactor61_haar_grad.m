function [f, gradf] = myfun_reactor61_haar_grad(w)
% Funcao objetivo COM GRADIENTE ANALITICO - Exemplo 6.1, versao Haar
% w = [tf; xdA(1:m); xdB(1:m); u(1:m)]
%
% min tf  -- trivial: gradiente e' [1; 0; 0; ...; 0]

    n = length(w);
    f = w(1);

    if nargout > 1
        gradf = zeros(n, 1);
        gradf(1) = 1.0;
    end

end
