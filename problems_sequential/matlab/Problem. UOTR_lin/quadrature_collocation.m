function [T, Y] = quadrature_collocation(odefun, tspan, y0, h, uz, ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal, est, uref, tref, nref)

          
% t,x,cell2mat(u)',(ks),ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,...
    %                     tal,est,optSPAIDER.prob.uref,...
    %                  optSPAIDER.prob.tref,optSPAIDER.prob.nref
    % tspan: intervalo de tempo [t0, tf]
    % y0: condição inicial (vetor de tamanho N)
    % h: tamanho do passo

    t0 = 0;
    tf = tspan(end);
    N = ceil((tf - t0) / h);  % Número de passos
    T = linspace(t0, tf, N + 1);  % Vetor de tempos
    M = length(y0);  % Número de EDOs
    Y = zeros(length(T), M);  % Matriz de soluções
    Y(1, :) = y0;  % Condição inicial

    % Pontos e pesos de Gauss-Legendre para 2 pontos
    c = [0.5 - 0.5/sqrt(3); 0.5 + 0.5/sqrt(3)];  % Pontos
    w = [1; 1];  % Pesos

    for n = 1:N
        tn = T(n);
        yn = Y(n, :)';  % Transpor para vetor coluna

        % Estimativas de quadratura
        k = zeros(M, 2);
        for i = 1:2
            k(:, i) = modelo_eqX(tn + c(i) * h, yn + h * sum(w(1:i)' .* k(:, 1:i), 2),uz, ks,1); 

            % t,x,cell2mat(u)',(ks),ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,...
                    %     tal,est,optSPAIDER.prob.uref,...
                    %  optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks):(ti(ks+1)-ti(ks))/d_ode: ti(ks+1)], cell2mat(pzo), optODE


             % ( t_out, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg)                                                    
        end
        
        % Atualizar solução usando quadratura
        Y(n + 1, :) = (yn + h * sum(w' .* k, 2))';  % Transpor para vetor linha
    end
end