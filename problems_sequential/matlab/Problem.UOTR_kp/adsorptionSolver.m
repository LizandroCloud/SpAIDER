clear,clc
%% MAIN
% Parameters
L   = 0.02;                                                                 % Column length
t0  = 0;                                                                    % Initial Time
tf  = 20;                                                                   % Final time
dt  = 0.05;                                                                 % Time step
t   = t0:dt:tf;                                                             % Time vector
dz  = 5e-5;                                                                 % Mesh size
z   = 0:dz:L;                                                               % Mesh generation
n   = numel(z);
cF_1  = [0.1 1.5 5 10 50 100 250 500 1000];                                   % Different feed concentration values
cF_2  = 0.1*ones(size(cF_1));
cF    = [cF_1;cF_2];
% Solution and plotting
for k = 1:size(cF,2)
    sol(k)  = adsorptionSolver(t,z,cF(:,k));
    name    = sprintf('Feed concentration of compound #1 cFeed = %f.\n',cF(1,k));
    figure
    subplot(4,2,1)
    plot(z,sol(k).y(1:n,1:5:end))
    xlabel('z')
    ylabel('c1')
    title('c1(z) at different t')
    subplot(4,2,2)
    plot(z,sol(k).y(n+1:2*n,1:5:end))
    xlabel('z')
    ylabel('q1')
    title('q1(z) at different t')
    subplot(4,2,3)
    plot(sol(k).x,sol(k).y(1:40:n,:))
    xlabel('t')
    ylabel('c1')
    title('c1(t) at different z')
    subplot(4,2,4)
    plot(sol(k).x,sol(k).y(n+1:40:2*n,:))
    xlabel('t')
    ylabel('q1')
    title('q1(t) at different z')
    subplot(4,2,5)
    plot(z,sol(k).y(2*n+1:3*n,1:5:end))
    xlabel('z')
    ylabel('c2')
    title('c2(z) at different t')
    subplot(4,2,6)
    plot(z,sol(k).y(3*n+1:4*n,1:5:end))
    xlabel('z')
    ylabel('q2')
    title('q2(z) at different t')
    subplot(4,2,7)
    plot(sol(k).x,sol(k).y(2*n+1:40:3*n,:))
    xlabel('t')
    ylabel('c2')
    title('c2(t) at different z')
    subplot(4,2,8)
    plot(sol(k).x,sol(k).y(3*n+1:40:4*n,:))
    xlabel('t')
    ylabel('q2')
    title('q2(t) at different z')
    sgtitle(name)
    
    disp('-------------------------------------------------------------------------------------------------')
end

%% Solver function
function out = adsorptionSolver(t,z,cFeed)
n = numel(z);                                                               % Size of mesh grid
h = diff(z); h = h(1);
% Initial Conditions
c1_0      = zeros(n,1);     
c1_0(1)   = cFeed(1);
q1_0      = zeros(n,1);
c2_0      = zeros(n,1);     
c2_0(1)   = cFeed(2);
q2_0      = zeros(n,1);
y0      = [c1_0 ; q1_0; c2_0 ; q2_0];
% Creating mass matrix
M               = eye(4*n);                                                 % Mass matrix
M(1,1)          = 0;                                                        % Algebraic equation for c1 at z = 0;
M(n,n)          = 0;                                                        % Algebraic equation for c1 at z = L
M(2*n+1,2*n+1)  = 0;                                                        % Algebraic equation for c2 at z = 0;
M(3*n,3*n)      = 0;                                                        % Algebraic equation for c2 at z = L
opts = odeset('Mass',M,'MassSingular','yes');
% Solving
out = ode15s(@(t,y) adsorptionModel(t,y,h,n,cFeed),t,y0,opts);
 
end
%% Adsorption model
function dydt = adsorptionModel(~,y,h,n,cFeed)
% Constants
D   = [1.85e-5;1.85e-5];                                                    % Axial Dispersion coefficient
v   = 3.38e-3 ;                                                             % Superficial velocity
epl = 0.47;                                                                 % Void fraction
k   = [1.9e-6;1.9e-6];                                                               % Mass Transfer Coefficient
b   = [0.05;0.1];                                                           % Langmuir parameter
qm  = [14.5;14.5];                                                          % Saturation capacity
% Variables allocation
c1 = y(1:n);
q1 = y(n+1:2*n);
c2 = y(2*n+1:3*n);
q2 = y(3*n+1:4*n);
% Langmuir term
c   = [c1,c2];
q   = [q1,q2];
LT  = (b'.*c)./(1+c*b)-q;
% Concentration equations
dc1dz            = zeros(n,1);
dc1dz(2:n-1)     = (c1(3:n)-c1(2:n-1))/h;
d2c1dz2          = zeros(n,1);
d2c1dz2(2:n-1)   = (c1(3:n)-2*c1(2:n-1)+c1(1:n-2))/h^2;
dc1dt            = -v*dc1dz+D(1)*d2c1dz2-k(1)*qm(1)*(1/epl-1)*LT(:,1);
dc1dt(1)         = c1(1)-cFeed(1);
dc1dt(n)         = c1(n)-c1(n-1);
dc2dz            = zeros(n,1);
dc2dz(2:n-1)     = (c2(3:n)-c2(2:n-1))/h;
d2c2dz2          = zeros(n,1);
d2c2dz2(2:n-1)   = (c2(3:n)-2*c2(2:n-1)+c2(1:n-2))/h^2;
dc2dt            = -v*dc2dz+D(2)*d2c2dz2-k(1)*qm(1)*(1/epl-1)*LT(:,2);
dc2dt(1)         = c2(1)-cFeed(2);
dc2dt(n)         = c2(n)-c2(n-1);
% Saturation equations
dq1dt = k(1)*qm(1)*LT(:,1);
dq2dt = k(2)*qm(2)*LT(:,2);
% Concatenate vectors of time derivatives
dydt = [dc1dt;dq1dt;dc2dt;dq2dt];
end