% Definir intervalo de tempo e condição inicial
tf=10;
tspan = [0, tf];
x0=[0.1;1;0.3]';
y0 = x0;  % Condição inicial para 4 EDOs

tic
I=[];
O=[];
for j = 4:20
u0=j; tx=[]; yx=[]; tend=0; ux0=[];
for k=1:10
    u=u0+2*(-1)^(k);
    [t, y] = ode15s(@systemODE, tspan/10, y0, [], u);
    tx=[tx;t(1:end)+tend];
    yx=[yx;y(1:end,:)];
    ux=[ux0; ones(length(t),1)*u];
    ux0=ux;
    y0=y(end,:);
    tend = tx(end);
end
I = [I; ux ];
O = [O; yx ];
end
toc
figure;
hold on;
for i = 1:size(yx, 2)
    plot(tx, yx(:, i), 'DisplayName', ['y' num2str(i)]);
end

xlabel('Tempo');
ylabel('Solução');
legend;
title('Solução do Sistema de 4 EDOs usando ode45');

% Criar conjunto de dados
data = iddata(O, I, 0.1); % 0.1 é o intervalo de amostragem

% Criar o modelo NLARX
na = 2; % Ordem da saída
nb = 2; % Ordem da entrada
nk = 0; % Atraso
orders = [3 1 1 1 0; % Ordem da saída
          2 1 1 1 0;
          2 1 1 1 0]; % Ordem da entrada
model = nlarx(data, orders);

% Validação do modelo
compare(data, model)
y_simulated = sim(model, data);

u0=15; tx=[]; yx=[]; tend=0; uz0=[];
x0=[0.1;1;0.3]';
y0 = x0;  % Condição inicial para 4 EDOs
tu = linspace(0,tf,10);
for k=1:9
    u=u0+2*(-1)^(k);
    % u = u0 - k*0.5;
    [t, y] = ode15s(@systemODE, [tu(k):(tu(k+1)-tu(k))/20:tu(k+1)], y0, [], u);
    tx=[tx;t(1:end)+tend];
    yx=[yx;y(1:end,:)];
    uz=[uz0; ones(length(t),1)*u];
    uz0=uz;
    y0=y(end,:);
    tend = tx(end);

end
nt = linspace(0,tf,length(uz));
data = iddata(yx, uz, 0.1); % 0.1 é o intervalo de amostragem
y_simulated = sim(model, data);
figure(10);
plot(nt,y_simulated.y(:,1))
hold on;
plot(nt,yx(:,1),'--')



function dx = systemODE(t, x, u)
    dx(1) = -x(1)^2 + u;
    dx(2) = -x(2)/3 + x(1);
    dx(3) = -x(2)^0.5/5 + u/5;
    dx=dx';
end

