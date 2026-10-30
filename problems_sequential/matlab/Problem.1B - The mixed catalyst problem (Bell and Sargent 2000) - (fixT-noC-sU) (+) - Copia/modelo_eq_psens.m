function [dx] = modelo_eq_psens(t_out,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est, yi,utarg,tstarg,nstarg)

global modelo_u0 modelo_par
global idu


for j=1:n_c
    ws_reg(j)=ver_reg_est1(ws(j), yi,t_out);
    % ws_reg(j) = get_u_smooth(t_out, te(j,2:ns_a(j)), ws(j,:) );
    % ws_reg(j) = u_reg1(t_out,te(j,:),ws(j),ks,ns_a(j));
    % ws_reg=ver_reg_est1(ws(j),yi,t_out);
    % ws_reg(j) = interp1(te(j,:), [ws(j,:),ws(j,end)], t_out, 'linear');
    % wsb_reg(j) = u_reg1(t_out,te(j,:),ws(j,:),ks,ne(j));
    % det=ws_reg(j) - wsb_reg(j);
    % if det>0.01
    %     det
    % end
for j=1:n_c
    u(j)=ws_reg(j);
end
end
dx(1)=(u(1)*(10*x(2)-x(1))              )*tal;
dx(2)=(u(1)*(x(1)-10*x(2))-(1-u(1))*x(2))*tal;
dx=transp(dx);

end