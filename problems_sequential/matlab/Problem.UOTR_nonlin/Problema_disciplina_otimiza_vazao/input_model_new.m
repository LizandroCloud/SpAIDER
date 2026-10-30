function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model


[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 1000; % Final time for integration
input.x0 = [ 3.24 1.44 0.1 0.01 0.001 0.001 0.888 0.036 0 0.006 0.0003 0.0063 0.0063 0.0063 0.0113 0.5]; %state Initial Condition
input.u_lb = [0   0   0   0]; % control vector lower boundary 
input.u_ub = [1.6 0.1 1.8 2]; % control vector upper boundary 
input.u0 = [0.08 0.001 0.09 0.001]; % scalar initial estimative of control profile for each variable
input.nU = 4; % Number of control variables...
input.is0 = 0;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 0;  % Final Iteration 


input.x0_lb = [0 0 0 ]; %x0 lb
input.x0_ub = [0.74 0.07 1.22]; %x0 lb

%--------------------------------------------------------------------------

% parameters
par(1)={'r1amax=0.85'}; 
par(2)={'r1bmax=0.022'};
par(3)={'r2mmax=0.018'};
par(4)={'r3max=0.09'};
par(5)={'r4max=0.03'};
par(6)={'r5max=0.03'};
par(7)={'r6max=0.75'};
par(8)={'r7max=0.03'};
par(9)={'r7mmax=0.012'};
par(10)={'k1a=0.08'};
par(11)={'k1b=0.001'};
par(12)={'k2m=0.009'};
par(13)={'k3=0.09'};
par(14)={'k4=0.05'};
par(15)={'k5=0.05'};
par(16)={'k6=0.05'};
par(17)={'k7=0.001'};
par(18)={'k7m=0.03'};
par(19)={'y1=1'};
par(20)={'y2m=1'};
par(21)={'y3=0.68'};
par(22)={'y4=0.7'};
par(23)={'y5=1'};
par(24)={'y6=0.7'};
par(25)={'y7=1.1'};
par(26)={'ke1a=20*k1a'};
par(27)={'remax=0.001'};
par(28)={'beta=0.1'};
par(29)={'rgmaxchapeu=0.0575'};
par(30)={'kglinha=0.01'};
par(31)={'kg=0.05'};
par(32)={'D = 0.031'};
par(33)={'ycm1=0.2'};
par(34)={'ycm2=0.15'};
par(35)={'ycm3=0.15'};
par(36)={'ycm4=(1-ycm1-ycm2-ycm3)/2'}; 
par(37)={'s1f=0.95'};
par(38)={'s2f=0.4238'};
par(39)={'alfa=100'};


%alg equations
y(1) = {'r1amax*(x(1)/(k1a+x(1)+((x(1)^2)/u(1))))*(x(7)/(x(7)+kglinha))*(x(10)/(((10^-7)+remax)/(rgmaxchapeu+beta)))'}; %r1a
y(2) = {'r1bmax*(x(1)/(k1b+x(1)+((x(1)^2)/u(2))))*(x(7)/(x(7)+kglinha))*(x(11)/(((10^-7)+remax)/(rgmaxchapeu+beta)))'};%r1b
y(3) = {'r3max*(x(2)/(k3+x(2)+((x(2)^2)/u(3))))*(x(7)/(x(7)+kglinha))*(x(12)/(((10^-7)+remax)/(rgmaxchapeu+beta)))'};%r3
y(4) = {'r7max*(x(2)/(k7+x(2)+((x(2)^2)/u(4))))*(x(7)/(x(7)+kglinha))*(x(15)/(((10^-7)+remax)/(rgmaxchapeu+beta)))'};%r7
y(5) = {'r4max*(x(3)/(k4+x(3)))*(x(7)/(x(7)+kglinha))*(x(13)/(((10^-7)+remax)/(rgmaxchapeu+beta)))'};%r4
y(6) = {'r6max*(x(3)/(k6+x(3)))*(x(7)/(x(7)+kglinha))*(x(14)/(((10^-7)+remax)/(rgmaxchapeu+beta)))'};%r6
y(7) = {'r5max*(x(4)/(k5+x(4)))*(x(7)/(x(7)+kglinha))'};%r5
y(8) = {'r2mmax*(x(1)/(k2m+x(1)))*(x(7)/(x(7)+kglinha))'};%r2m
y(9) = {'r7mmax*(x(2)/(k7m+x(2)))*(x(7)/(x(7)+kglinha))'};%r7m
y(10) = {'remax*(x(1)/(ke1a+x(1)))'};%re1a
y(11) = {'remax*(x(1)/(k1b+x(1)))'};%re1b
y(12) = {'remax*(x(2)/(k3+x(2)))'};%re3
y(13) = {'remax*(x(3)/(k4+x(3)))'};%re4
y(14) = {'remax*(x(3)/(k6+x(3)))'};%re6
y(15) = {'remax*(x(2)/(k7+x(2)))'};%re7
y(16) = {'rgmaxchapeu*(x(3)/(x(3)+kg))*(x(4)/(x(4)+kg))*(x(5)/(x(5)+kg))*(x(6)/(x(6)+kg))'};%rgchapeure
y(17)={'(1+tanh((y(1)-y(2))*alfa))/2'};
y(18)={'(1+tanh((y(3)-y(4))*alfa))/2'};
y(19)={'(1+tanh((y3*y(3)-y4*y(5))*alfa))/2'};
y(20)={'y(1)/(y(1)*y(17)+y(2)*(1-y(17)))'}; 
y(21)={'y(2)/(y(1)*y(17)+y(2)*(1-y(17))'}; 
y(22)={'y(3)/(y(3)*y(18)+y(4)*(1-y(18)))'};
y(23)={'y(4)/(y(3)*y(18)+y(4)*(1-y(18)))'};
y(24)={'y3*y(3)/((y3*y(3))*y(19)+(y4*y(5))*(1-y(19)))'}; 
y(25)={'y4*y(5)/((y3*y(3))*y(19)+(y4*y(5))*(1-y(19)))'};
y(26)={'y1*(y(1)*y(20)+y(2)*y(21))+y2m*y(8)+y3*y(3)*y(22)+y7*(y(4)*y(23)+y(9))+(y4-1)*y(5)*y(25)+(y5-1)*y(7)-y(6)*y(20)+y(16)*(1-ycm1-ycm2-ycm3-ycm4)'};
y(27)={'y(1)/(y(1)+y(2))'}; 
y(28)={'y(2)/(y(1)+y(2))'}; 
y(29)={'y(3)/(y(3)+y(4))'};
y(30)={'y(4)/(y(3)+y(4))'}; 
y(31)={'y3*y(3)/(y3*y(3)+y4*y(5))'};
y(32)={'y4*y(5)/(y3*y(3)+y4*y(5))'};   


% ode equations
dx(1) = {'-(y(1)*y(20) + y(2)*y(21) +y(8))*x(8) +D*(s1f-x(1))'};%ds1/dt
dx(2) = {'-(y(3)*y(24)*y(22)+y(4)*y(23)+y(9))*x(8)+D*(s2f-x(2)) '};%ds2/dt
dx(3) = {'y1*(y(1)*y(20)+y(2)*y(21))+y5*y(7)-y(5)*y(25)-y(6)*y(20)-ycm1*y(16)-y(26)*x(3)'};%dm1/dt
dx(4) = {'y3*y(3)*y(24)*y(22)+y4*y(5)*y(25)-y(7)-ycm2*y(16)-y(26)*x(4)'};%dm2/dt
dx(5) = {'y2m*y(8)-ycm3*y(16)-y(26)*x(5)'};%dm3/dt
dx(6) = {'y7*(y(4)*y(23)+y(9))-ycm4*y(16)-y(26)*x(6)'};%dm4/dt
dx(7) = {'y(16)-y(26)*x(7)'};%dc'/dt
dx(8) = {'(y(26)-D)*x(8)'};%dx/dt
dx(9) = {'y6*y(6)*y(20)*x(8)-D*x(9)'};%dL/dt
dx(10) = {'(10^-7)+y(10)*y(27)-(y(26)+beta)*x(10)'};%de1a/dt
dx(11) = {'(10^-7)+y(11)*y(28)-(y(26)+beta)*x(11)'};%de1b/dt
dx(12) = {'(10^-7)+y(12)*y(31)*y(29)-(y(26)+beta)*x(12)'};%de3/dt
dx(13) = {'(10^-7)+y(13)*y(32)-(y(26)+beta)*x(13)'};%de4/dt
dx(14) = {'(10^-7)+y(14)*y(27)-(y(26)+beta)*x(14)'};%de6/dt
dx(15) = {'(10^-7)+y(15)*y(30)-(y(26)+beta)*x(15)'}; %de7/dt
dx(16) = {'0'};%dV/dt


% dx(4) = {'(u(1) - uref)^2'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'(-f(t_out,1) - 1) '};   % calculo da Fobj no tempo final...

% constraints

cr_eq(1) = {'[]'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
%cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality
%cr_deq(2) = {'x(end,2) - cb_mais'} ; % inequality

end