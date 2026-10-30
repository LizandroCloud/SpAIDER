function [dx] = modelo_eq_psens(t_out,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est, yi,utarg,tstarg,nstarg)

global modelo_u0 modelo_par
global idu

% idu=[1;16;17;18];
% idu=[1;18];
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

end
%Input interpolation
uf=modelo_u0';
for j=1:n_c
    uf(idu(j))=ws_reg(j);
end
uf(idu(1))=uf(idu(1))*1e3;

%States' derivative
dx=uotr_ctrl_Kp(t_out,x(1:n_s),uf,1,[],modelo_par);
dx=real(dx)*tal;

end