function J=minimiza(par)

global x0 M1 M2 modelo t_i t_f c_m contador

contador=contador+1;

try

op=odeset('MaxStep',0.01);
op1=odeset('AbsTol',1e-5);
[t,x]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,par),[0:1:t_i],x0);
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,par),[t_i:0.01:t_f],x(end,:),op);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,par),[t_f 11],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,par),[11:0.00001:15],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,par),[15:0.00001:18],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,par),[18:0.00001:20],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,par),[20:0.00001:21],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,par),[21:0.00001:22],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];
 [ti,xi]=ode15s(@(t,x) pseudo_steady_state(t,x,modelo,par),[22:0.00001:23],x(end,:),op1);
 x=[x(1:end-1,:);xi]; t=[t(1:end-1);ti];

catch ME
  J = 1e30;
  return
end

Enz=pseudo_steady_state_v(x);
p_o2=100*x(:,9)/c_m;

x_1=[x(:,2) x(:,3) x(:,6) x(:,4) p_o2 Enz x(:,5)];

g=[];
for i=1:7
g(:,i) = x_1(:,i); 
t1{1,i} = M1{1,i}; 
x1{1,i} = interp1(t,g(:,i),t1{1,i},'spline'); 
end

for j=1:7
    S=0;
        for i=1:length(x1)
            S = S+(x1{1,j}(i) - M2{1,j}(i))^2;
        end
    Z(j) = S;
end

J=sum(Z);
save([num2str(contador),'.mat'],'J','par')

