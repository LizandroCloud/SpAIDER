function [fobj] = openloop(L)
N=300; % 3*(# of spatial discretizations)
x=linspace(0,1,N/3); % Discretizing the spatial axis for plotting
uinitial=zeros(N,1); % Initialing the Initial condition vector
tspan = [0 7000]; 
% Initial condition for dimensionless concentration for the reaction stream
for i=1:N/3
uinitial(i)= 1; %was 1
end
% Initial condition for dimensionless temp for the reaction stream
for i=N/3+1:2*N/3
uinitial(i)=0;
end
% Initial condition for dimensionless temperature for the wall
for i=(2*N)/3+1:N
uinitial(i)=0;
end
uinitial = uinitial';
% Using ode solver to solve the system of discretized pde
[T,Y] = ode23s(@(t,u) sub(t,u,L),tspan,uinitial);
[a, b]=size(Y);
% Converting the output matrix from ode solver to give us only dimensionless concentration with time and spatially
u=zeros(a,N/3);
for i=1:a
 for j = 1:(b/3)
 u(i,j)= Y(i,j);
 end 
end
save('heatingcoilconcinitial.mat','u') % Saving the matrix for IC of bvp4c code
% Converting the output matrix from ode solver to give us only dimensionless temperature with time and spatially
v=zeros(a,N/3);
for i=1:a
 for j = (b/3+1):(2*b)/3
 v(i,j-b/3)= Y(i,j);
 end
end
save('heatingcoiltempinitial.mat','v') % Saving the matrix for IC of bvp4c code
% Converting the output matrix from ode solver to give us only dimensionless temperature for the wall with time and spatially
vw=zeros(a,N/3);
for i=1:a
 for j = ((2*b)/3+1):b
 vw(i,j-(2*b)/3)= Y(i,j); 
 end
end
save('heatingcoiltempwinitial.mat','vw') % Saving the matrix for IC bvp4c code
% plot(T,1-u(:,N/3))
Outlet_stream_conv = 1-u(a,:);
Outlet_stream_temp = v(a,N/3);
Inlet_wall_temp = vw(a,2);
fobj = u(end,100);
end