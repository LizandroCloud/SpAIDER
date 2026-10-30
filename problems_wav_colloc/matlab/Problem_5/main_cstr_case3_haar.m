% Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
% Case III (Balsa-Canto et al., 2001) - "Non-linear CSTR", versao Haar (Hi/iDWT)
%
% Mesma fisica do main_cstr_case3.m, mas reconstruindo TODOS os 8
% estados via a matriz de Haar (H_M.m/P_M.m/wav_mat.m), a partir de
% suas derivadas amostradas -- ao contrario do dominio fisico, aqui
% NENHUM estado e' uma variavel de decisao direta.
%
% w = [xd1(1:m), ..., xd8(1:m), u1(1:m), ..., u_nu(1:m)]
%
% Requer no path: H_M.m, P_M.m, wav_mat.m, haar.m, haar2.m
%
% *** NOTA IMPORTANTE (ver myfun_cstr_case3_haar.m) ***
% Esse sistema tem controles otimos de forma bastante "selvagem" (nas
% palavras do proprio artigo). O operador Hi/iDWT (global) converge
% mais devagar que a colocacao trapezoidal do dominio fisico para esse
% tipo de dinamica -- validado numericamente (Python): m=16 -> J*~19.57,
% m=32 -> J*~19.83, convergindo rumo a referencia 20.09 mas precisando
% de mais refinamento que main_cstr_case3.m para a mesma precisao. Se
% voce so precisa do resultado numerico (nao da base de Haar em si),
% prefira main_cstr_case3.m.
%
% Referencias de otimo:
%   3 controles: J* = 20.0895 (rho=10), 20.0953 (rho=40)
%   4 controles: J* = 21.757  (rho=11), 21.807 (rho=80)

clear; close all; clc;

z0 = [0.1883, 0.2507, 0.0467, 0.0899, 0.1804, 0.1394, 0.1046, 0.000];
tf = 0.2;
u_lb = [0.0, 0.0, 0.0, 0.0];
u_ub = [20.0, 6.0, 4.0, 20.0];
u4_fixed = 6.0;

lambda_reg = 0.001;

optNLP = optimset('Algorithm', 'sqp', 'Display', 'off', ...
    'TolX', 1e-11, 'TolFun', 1e-11, 'TolCon', 1e-11, ...
    'MaxFunEval', 300000, 'MaxIter', 6000);

%**************************************************
% Sub-caso (i): 3 controles (u4 fixo em 6.0)
%**************************************************
fprintf('=== Sub-caso 3 controles (Haar) ===\n');
n_u = 3;

J = 3;   % m=16
m = 2^(J+1);
[H, ~, Hi] = H_M(J);
[~, ~, PH] = P_M(J, H);
[~, iDWT]  = wav_mat(PH, Hi, J);

u0 = reshape(repmat((u_lb(1:n_u)+u_ub(1:n_u))'/2, 1, m)', 1, []);
w0 = [zeros(1, 8*m), u0];
lb = [-Inf(1, 8*m), reshape(repmat(u_lb(1:n_u)', 1, m)', 1, [])];
ub = [ Inf(1, 8*m), reshape(repmat(u_ub(1:n_u)', 1, m)', 1, [])];

[w, fval, exitflag] = fmincon(@(w) myfun_cstr_case3_haar(w, m, n_u, Hi, iDWT, z0, lambda_reg), w0, ...
    [], [], [], [], lb, ub, @(w) mycon_cstr_case3_haar(w, m, n_u, Hi, iDWT, z0, tf, u4_fixed), optNLP);

A = Hi \ iDWT;
xd8 = w(7*m+1:8*m);
z8 = z0(8) + xd8*A;
fprintf('J=%d (m=%d): J* = %.4f  | exitflag=%d\n', J, m, z8(end), exitflag);

w_prev = w; m_prev = m;
for J2 = [4:1 :6]
    m2 = 2^(J2+1);
    rep = m2/m_prev;

    w0 = [];
    for i = 1:8
        xdi_prev = w_prev((i-1)*m_prev+1 : i*m_prev);
        w0 = [w0, repelem(xdi_prev, rep)]; %#ok<AGROW>
    end
    for i = 1:n_u
        ui_prev = w_prev(8*m_prev + (i-1)*m_prev+1 : 8*m_prev + i*m_prev);
        w0 = [w0, repelem(ui_prev, rep)]; %#ok<AGROW>
    end

    [H2, ~, Hi2] = H_M(J2);
    [~, ~, PH2]  = P_M(J2, H2);
    [~, iDWT2]   = wav_mat(PH2, Hi2, J2);

    lb2 = [-Inf(1, 8*m2), reshape(repmat(u_lb(1:n_u)', 1, m2)', 1, [])];
    ub2 = [ Inf(1, 8*m2), reshape(repmat(u_ub(1:n_u)', 1, m2)', 1, [])];

    [w, fval, exitflag] = fmincon(@(w) myfun_cstr_case3_haar(w, m2, n_u, Hi2, iDWT2, z0, lambda_reg), w0, ...
        [], [], [], [], lb2, ub2, @(w) mycon_cstr_case3_haar(w, m2, n_u, Hi2, iDWT2, z0, tf, u4_fixed), optNLP);

    A2 = Hi2 \ iDWT2;
    xd8_2 = w(7*m2+1:8*m2);
    z8_2 = z0(8) + xd8_2*A2;
    [c, ceq] = mycon_cstr_case3_haar(w, m2, n_u, Hi2, iDWT2, z0, tf, u4_fixed);
    fprintf('J=%d (m=%d): J* = %.4f | exitflag=%d | max|ceq|=%.2e\n', J2, m2, z8_2(end), exitflag, max(abs(ceq)));

    w_prev = w; m_prev = m2;
end
fprintf('Referencia: J* = 20.0895 (rho=10), 20.0953 (rho=40)\n\n');

%**************************************************
% Plot (3 controles)
%**************************************************
m = m_prev;
A = Hi2 \ iDWT2;
U = zeros(n_u, m);
for i = 1:n_u
    U(i,:) = w(8*m + (i-1)*m+1 : 8*m+i*m);
end
t_nodes = ((1:m) - 0.5)/m*tf;

figure;
stairs(t_nodes, U(1,:), 'LineWidth', 1.2); hold on;
stairs(t_nodes, U(2,:), 'LineWidth', 1.2);
stairs(t_nodes, U(3,:), 'LineWidth', 1.2);
xlabel('t'); ylabel('u'); legend('u_1','u_2','u_3','Location','best');
title('Case III (Haar) - 3 controles');
