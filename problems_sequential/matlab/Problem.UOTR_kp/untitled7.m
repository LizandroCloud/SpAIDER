% Definir o intervalo de tempo e a condição inicial
tspan = [0, 0.01];
y0 = [1; 1];  % Condição inicial [y10, y20]
h = 0.0001;

% Chamar a função do método dos trapézios
[T, Y] = trapezoidal_method_2d(@systemODE, tspan, y0, h);

% Plotar os resultados
plot(T, Y(:,1), 'r', T, Y(:,2), 'b');
xlabel('Tempo');
ylabel('Solução');
legend('y1', 'y2');
title('Solução do Sistema de EDOs usando o Método dos Trapézios');



function [T, Y] = trapezoidal_method_2d(odefun, tspan, y0, h)
    % tspan: intervalo de tempo [t0, tf]
    % y0: condição inicial [y10, y20]
    % h: tamanho do passo

    t0 = tspan(1);
    tf = tspan(2);
    N = (tf - t0) / h;  % Número de passos
    
    T = t0:h:tf;
    Y = zeros(length(T), 2);  % Para duas variáveis
    Y(1, :) = y0;

     for n = 1:N
        tn = T(n);
        yn = Y(n, :)';
        tn1 = T(n+1);

        % Resolver as equações implícitas dos trapézios
        % Inicializar a solução para o próximo passo
        y_next = yn;  % Valor inicial para a iteração

        % Iteração de Newton para resolver o sistema implícito
        for iter = 1:10
            % Calcular os resíduos
            f_tn = odefun(tn, yn);
            f_tn1 = odefun(tn1, y_next);

            F1 = yn(1) + (h/2) * (f_tn(1) + f_tn1(1)) - y_next(1);
            F2 = yn(2) + (h/2) * (f_tn(2) + f_tn1(2)) - y_next(2);
            F = [F1; F2];
            
            % Jacobiana do sistema
            J = eye(2) - (h/2) * [1 -1; 1 1];
            y_next = y_next - J \ F;  % Iteração de Newton
        end
        
        % Armazenar a solução
        Y(n+1, :) = y_next';
    end
end



function dydt = systemODE(t, y)
    dydt = [0.5*y(2);  % dy1/dt
            -0.9*y(1)]; % dy2/dt
end

   % % Função de EDO fornecida
   %  function dydt = odefun(t, y)
   %      dydt = [y(1) - y(2);  % dy1/dt
   %              y(1) + y(2)]; % dy2/dt
   %  end

