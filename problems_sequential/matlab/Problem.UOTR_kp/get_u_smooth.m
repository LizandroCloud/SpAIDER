function u = get_u_smooth(t, ti, ui)
    % t    : tempo atual
    % ti   : vetor de tempos de troca
    % ui   : vetor de valores desejados
    % k    : parâmetro de suavidade (maior k = transições mais bruscas)
    k = 100;
    n = length(ti);
    u = 0;

    for i = 1:n-1
        s = 0.5 * (1 + tanh(k * (t - ti(i))));
        u = u + (ui(i+1) - ui(i)) * s;
    end

    u = u + ui(1);  % valor base
end
