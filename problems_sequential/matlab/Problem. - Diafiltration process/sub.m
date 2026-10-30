function dudt=sub(t,u,L)
N=300; % Definning the number of points in spatial axis
x=linspace(0,3,N); % Discretizing the spatial axis
deltax=x(2)-x(1); % Definning each interval size
dudt = zeros(N,1); % Initializing the vector for all 3 PDEs combined
% The boundary conditions are--
u(1)=1; u(N/3+1)= 0; % Inlet conditions for concentration and temperature
u((2*N)/3+1)= u((2*N)/3+2); u(N)=u(N-1); % Adiabatic boundary conditions for the wall
% The dimensionless groups are--
Da = 20.6; alpha= L; gamma= -0.636 ;NTU=188; CP=11.7; theta=0.25; 
Fo=0.00293; psi=-0.045;
% psi = 0.045+0.155*heaviside(10)+0.155*heaviside(40);
% Order of the reactions
n=1;
for i=2:N/3
 % Dimensionless concentration for the reaction stream
 dudt(i) = -(u(i)-u(i-1))/deltax -Da*theta*exp((alpha*gamma*u(N/3+i-1))/(1+gamma*u(N/3+i-1)))*u(i-1)^n; 
 % Dimensionless temp for the reaction stream
 dudt(N/3+i) = -(u(N/3+i)-u(N/3+i-1))/deltax + ...
Da*theta*exp((alpha*gamma*u(N/3+i-1))/(1+gamma*u(N/3+i-1)))*u(i-1)^n -...
NTU*theta*(u(N/3+i-1)-u((2*N)/3+i-1));
end
for i=2:N/3-1
% Dimensionless temperature for the wall
 dudt((2*N)/3+i) = Fo*theta*((u((2*N)/3+i+1)+u((2*N)/3+i-1)-...
2*u((2*N)/3+i))/deltax^2) + Fo*theta*(NTU/CP*(u(N/3+i-1)-u((2*N)/3+i-1))+psi);
end
end