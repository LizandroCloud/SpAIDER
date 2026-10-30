function f = myfun_cstr_case3_haar(w, m, n_u, Hi, iDWT, z0, lambda_reg)
% Funcao objetivo - Case III (CSTR nao linear), versao Haar (Hi/iDWT).
%
% w = [xd1(1:m), xd2(1:m), ..., xd8(1:m), u1(1:m), ..., u_nu(1:m)]
%   (cada bloco de tamanho m, concatenado -- mesma convencao dos
%   scripts anteriores, NAO interleado por tempo)
%   xd_i = amostras de dz_i/dtau (tempo normalizado tau=[0,1])
%   u_i  = controles amostrados nos mesmos pontos (n_u = 3 ou 4)
%
% tf e' FIXO (0.2). max z8(tf)  <=>  min f = -z8(tf) + regularizacao
%
% NOTA: esse sistema tem controles otimos de forma "bastante selvagem"
% (o proprio artigo de Balsa-Canto et al. menciona isso). O operador
% Hi/iDWT (global) converge mais devagar que a colocacao trapezoidal
% local para esse tipo de dinamica -- validado numericamente: m=16 dá
% J*~19.57, m=32 dá J*~19.83, convergindo rumo a referencia 20.09, mas
% precisando de malha mais fina que a versao em dominio fisico
% (main_cstr_case3.m) para a mesma precisao.

    if nargin < 7
        lambda_reg = 0.001;
    end

    A = Hi \ iDWT;

    xd8 = w(7*m+1 : 8*m);
    z8 = z0(8) + xd8*A;

    f = -z8(end);

    for i = 1:n_u
        ui = w(8*m + (i-1)*m + 1 : 8*m + i*m);
        f = f + lambda_reg*sum(diff(ui).^2);
    end

end
