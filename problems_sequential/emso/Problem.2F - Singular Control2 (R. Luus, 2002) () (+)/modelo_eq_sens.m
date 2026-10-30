function [ dx ] = modelo_eq_sens( t_out, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
[u] = ver_reg_est1(ws,yi,t_out);
[uref] = ver_reg_est1(ws,yi,t_out);
;


dx(1)=(x(2)  )*tal;
dx(2)=(u(1)  )*tal;
dx(3)=(x(1)^2)*tal;
dx=transp(dx);
jc1 = [
([0, 1, 0, 0])     
([0, 0, 0, 1])     
([2*x(1), 0, 0, 0])
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