function f = myfun_reactor61_haar(w, m, lambda_reg)
% Funcao objetivo - Exemplo 6.1, versao Haar (Hi/iDWT)
% w = [tf; xdA(1:m); xdB(1:m); u(1:m)]
%   tf  = tempo final (LIVRE)
%   xdA, xdB = amostras de cA_dot, cB_dot (em tempo normalizado tau)
%   u   = controle amostrado (vazao de B)
%
% min tf + lambda_reg * sum((u_k+1 - u_k)^2)
%
% O termo de regularizacao penaliza variacao brusca entre pontos
% vizinhos de u. Sem ele, o problema discretizado tem uma direcao
% "achatada" em u -- varios perfis diferentes de u dao o MESMO tf*
% (nao e' um arco singular do problema continuo, e' uma degenerescencia
% da discretizacao), e o otimizador pode parar num membro serrilhado
% dessa familia dependendo so do chute inicial de tf. lambda_reg
% pequeno (ex: 0.01) elimina o serrilhado SEM alterar tf* (validado
% numericamente).

    if nargin < 3
        lambda_reg = 0.01;
    end

    tf = w(1);
    u = w(2+2*m : 1+3*m);

    f = tf + lambda_reg * sum(diff(u).^2);

end
