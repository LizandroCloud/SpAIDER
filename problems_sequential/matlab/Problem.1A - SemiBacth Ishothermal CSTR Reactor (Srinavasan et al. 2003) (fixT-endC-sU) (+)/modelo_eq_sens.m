function [ dx ] = modelo_eq_sens( t_out, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
[u] = ver_reg_est1(ws,yi,t_out);
[uref] = ver_reg_est1(ws,yi,t_out);
k1 = 0.053     ;
k2 = 0.128     ;
cb_in = 5      ;
Vo = 1         ;
lamb=1^-8      ;
cao=0.72       ;
cbo=0.05       ;
cd_mais = 0.15 ;
cb_mais = 0.025;


dx(1)=(-k1*x(1)*x(2)-( u(1)/x(3) )*x(1)                           )*tal;
dx(2)=(-k1*x(1)*x(2) - 2*k2*x(2)^2 + ( (u(1)/x(3))*(cb_in - x(2))))*tal;
dx(3)=(u(1)                                                       )*tal;
dx=transp(dx);
jc1 = [
[- k1*x(2) - u(1)/x(3), -k1*x(1), (u(1)*x(1))/x(3)^2, -x(1)/x(3)]                                
[-k1*x(2), - k1*x(1) - 4*k2*x(2) - u(1)/x(3), -(u(1)*(cb_in - x(2)))/x(3)^2, (cb_in - x(2))/x(3)]
[0, 0, 0, 1]                                                                                     
]*tal;
jk=double(jc1);
j0=0;
jk(:,1:n_s)=jc1(:,1:n_s);
jk(:,n_s+1:n_s+n_c)=0;
for j=1:n_c
for is = 1:ns_a(j)
jk(:,n_s+j)=uws(j,is)*jc1(:,n_s+j);
k=1:(n_s);
dx(n_s*is+k+j0) = jk(k,1:n_s)*x(n_s*is+1+j0:n_s*is+n_s+j0)+ jk(k,n_s+j);
end 
j0=length(dx)-n_s;
end 
end