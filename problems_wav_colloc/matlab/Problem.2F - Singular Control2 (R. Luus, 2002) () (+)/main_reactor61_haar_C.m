% Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
% Exemplo 6.1 (Srinivasan, Palanki & Bonvin, 2003) - versao Haar (Hi/iDWT)
% "Isothermal semi-batch reactor with a safety constraint"
%
% Mesma fisica do main_reactor61.m (dominio fisico), mas agora
% reconstruindo cA(t), cB(t), V(t) via a matriz de Haar (H_M.m/P_M.m/
% wav_mat.m), no mesmo espirito do main_case_6_3_haar.m.
%
% w = [tf; xdA(1:m); xdB(1:m); u(1:m)]
%   tf       = tempo final LIVRE
%   xdA,xdB  = amostras de cA_dot, cB_dot (tempo normalizado tau=[0,1])
%   u        = vazao de B amostrada nos mesmos pontos de colocacao
%
% Requer no path: H_M.m, P_M.m, wav_mat.m, haar.m, haar2.m
%
% Referencia de otimo: tf* = 19.80 h, ts = 11.44 h (troca upath->umin)

clear; close all; clc;

%**************************************************
% Parametros (Tabela 3 do artigo)
%**************************************************
p.k    = 0.0482;
p.DH   = 60000.0;
p.rho  = 900.0;
p.cp   = 4.2;
p.cBin = 2.0;
p.umin = 0.0;
p.umax = 0.1;
p.Tmax = 80.0;
p.Vmax = 1.0;
p.nCdes = 0.6;
p.cAo  = 2.0;
p.Vo   = 0.7;
T      = 70.0;

p.cBmax = p.rho*p.cp*(p.Tmax - T) / p.DH;
p.cBo   = p.cBmax;

lambda_reg = 0.01;   % peso da regularizacao de suavidade em u (ver myfun_reactor61_haar.m)

fprintf('cBmax = %.4f mol/l\n', p.cBmax);

optNLP = optimset('Algorithm', 'sqp', 'Display', 'off', ...
    'TolX', 1e-10, 'TolFun', 1e-10, 'TolCon', 1e-10, ...
    'MaxFunEval', 200000, 'MaxIter', 5000);

%**************************************************
% Passo 1: malha grosseira (J=3 -> m=16), chute com tf perto do esperado
%**************************************************
J = 3;
m = 2^(J+1);

[H, ~, Hi] = H_M(J);
[~, ~, PH] = P_M(J, H);
[~, iDWT]  = wav_mat(PH, Hi, J);

w0 = [19.8, zeros(1,m), zeros(1,m), 0.05*ones(1,m)];
lb = [0.1,  -Inf(1,m), -Inf(1,m), p.umin*ones(1,m)];
ub = [100,   Inf(1,m),  Inf(1,m), p.umax*ones(1,m)];

[w, fval, exitflag, output] = fmincon(@(w) myfun_reactor61_haar_C(w, m, lambda_reg), w0, ...
    [], [], [], [], lb, ub, @(w) mycon_reactor61_haar(w, m, Hi, iDWT, p), optNLP);

fprintf('J=%d (m=%d): tf* = %.4f h  | exitflag=%d\n', J, m, w(1), exitflag);

%**************************************************
% Passo 2: refinamento de malha, reaproveitando a solucao anterior
%**************************************************
w_prev = w;
m_prev = m;
for J2 = [4 5]
    m2 = 2^(J2+1);
    rep = m2 / m_prev;

    tf_prev  = w_prev(1);
    xdA_prev = w_prev(2 : 1+m_prev);
    xdB_prev = w_prev(2+m_prev : 1+2*m_prev);
    u_prev   = w_prev(2+2*m_prev : 1+3*m_prev);

    w0 = [tf_prev, repelem(xdA_prev, rep), repelem(xdB_prev, rep), repelem(u_prev, rep)];

    [H2, ~, Hi2] = H_M(J2);
    [~, ~, PH2]  = P_M(J2, H2);
    [~, iDWT2]   = wav_mat(PH2, Hi2, J2);

    lb2 = [0.1,  -Inf(1,m2), -Inf(1,m2), p.umin*ones(1,m2)];
    ub2 = [100,   Inf(1,m2),  Inf(1,m2), p.umax*ones(1,m2)];

    [w, fval, exitflag, output] = fmincon(@(w) myfun_reactor61_haar_C(w, m2, lambda_reg), w0, ...
        [], [], [], [], lb2, ub2, @(w) mycon_reactor61_haar(w, m2, Hi2, iDWT2, p), optNLP);

    [c, ceq] = mycon_reactor61_haar(w, m2, Hi2, iDWT2, p);
    u_diag = w(2+2*m2 : 1+3*m2);
    fprintf('J=%d (m=%d): tf* = %.4f h  | exitflag=%d | max|ceq|=%.2e | max(c)=%.2e | max|du|=%.4f\n', ...
        J2, m2, w(1), exitflag, max(abs(ceq)), max(c), max(abs(diff(u_diag))));

    w_prev = w;
    m_prev = m2;
end

fprintf('\nReferencia (Srinivasan et al., 2003): tf* = 19.80 h,  ts = 11.44 h\n');

%**************************************************
% Reconstrucao final e plot
%**************************************************
m = m_prev;
[H, ~, Hi] = H_M(log2(m) - 1);
[~, ~, PH] = P_M(log2(m) - 1, H);
[~, iDWT]  = wav_mat(PH, Hi, log2(m) - 1);

tf  = w(1);
xdA = w(2 : 1+m);
xdB = w(2+m : 1+2*m);
u   = w(2+2*m : 1+3*m);

Hi_inv = inv(Hi); %#ok<MINV>
cA = p.cAo + (xdA * Hi_inv) * iDWT;
cB = p.cBo + (xdB * Hi_inv) * iDWT;
V  = p.Vo  + tf * ((u * Hi_inv) * iDWT);
cC = (p.cAo*p.Vo - cA.*V) ./ V;

t_nodes = ((1:m) - 0.5) / m * tf;   % pontos de colocacao, em tempo fisico

figure;
subplot(2,1,1);
stairs(t_nodes, u, 'LineWidth', 1.2);
ylabel('u [l/h]'); title(sprintf('Exemplo 6.1 (Haar) - vazao de B (tf=%.2f h)', tf));

subplot(2,1,2);
plot(t_nodes, cA, '-', t_nodes, cB, '-', t_nodes, cC, '-', 'LineWidth', 1.2);
hold on; yline(p.cBmax, '--k');
xlabel('t [h]'); ylabel('concentracao [mol/l]');
legend('c_A', 'c_B', 'c_C', 'c_{B,max}', 'Location', 'best');
