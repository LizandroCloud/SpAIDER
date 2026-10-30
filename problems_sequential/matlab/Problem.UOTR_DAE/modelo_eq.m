function dx = modelo_eq(t_out,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg)
global modelo_u0 modelo_par

for j=1:n_c
ws_reg=ver_reg_est0(ws,yi{j},t_out);
% ws_reg(j) = u_reg1(t_out,te(j,:),nstarg(:,j),ks,ns_a(j));
end

%idx=[2,3,5,6,7,9,10,13,15,17,19,20,22,24,25,28];
%idy=[1,4,8,12,14,16,18,26]; %idy=[1,4,8,11,12,14,16,18,21,23,26,27];

idu=[1]; %idu=[1,5,7,9,14,15,16,19,20,21,23,25,26];
%idw=[2,3,4,6,8,10,11,12,13,17,18,22,24,27];

u=modelo_u0; u(idu)=ws_reg*1e3;
dx=uotr_ctrl_dae(t_out,x,u,1,[],modelo_par);
dx=dx*tal;

end