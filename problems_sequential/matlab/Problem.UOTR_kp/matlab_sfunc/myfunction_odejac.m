%function [dFdx,dFdu] = myfunction_odejac(t_out,x,ws,ts,n_s,n_c,dx,dx_old,dFdx,dFdu)
function Jac = myfunction_odejac(t_out,x,ws,ts,n_s,n_c,ns_a)
global modelo_u0 modelo_par
global modelo_dFdx

%Nx=size(modelo_x0,1); Nu=size(modelo_u0(1),1);
%Nx=n_s; Nu=n_c;

uf=modelo_u0; uf(1)=interpl(t_out,ws,ts,1e-2)*1e3; uf=uf(:);

%Update dFdx
fhand=@(t,x)uotr_ctrl(t,x,uf,1,[],modelo_par);
fval=fhand([],x(1:n_s));
if imag(fval)>=1e-6
    error('Solution of states is imaginary')
end                    
[dFdx,~]=numjac(fhand,[],x(1:n_s),real(fval),zeros(n_s,1));
if imag(dFdx)>=1e-6
    error('Solution of states is imaginary')
end            
dFdx=real(dFdx);

%ODE system states + sensitivity jacobian
Jac=dFdx;
for j=1:n_c
    for is = 1:ns_a(j)
        %Jac=blkdiag(Jac,dFdx);
        Jac=blkdiag(Jac,modelo_dFdx);
    end
end

%Update dFdu
%fhand=@(t,u)uotr_ctrl(t,x,[u;modelo_u0(2:end)],1,[],modelo_par);
%fval=fhand([],uf(1));
%if imag(fval)>=1e-6
%    error('Solution of states is imaginary')
%end                            
%[dFdu,~]=numjac(fhand,[],uf(1),real(fval),zeros(n_c,1));
%dFdu=real(dFdu);

%{
for k=2:length(t_out)
    uf=modelo_u0; uf(1)=interpl(t_out(k),ws,ts,1e-2)*1e3; uf=uf(:);
    if abs((dx(k,:)'-dx_old(k,:)')./dx(k,:)')>1e-3    %for abrupt change in dx
        %Update dFdx
        fhand=@(t,x)uotr_ctrl(t,x,uf,1,[],modelo_par);
        fval=fhand([],x(k,:)');
        if imag(fval)>=1e-6
            error('Solution of states is imaginary')
        end                    
        [dFdx(:,:,k),~]=numjac(fhand,[],x(k,:)',real(fval),zeros(Nx,1));
        if imag(dFdx(:,:,k))>=1e-6
            error('Solution of states is imaginary')
        end            
        dFdx(:,:,k)=real(dFdx(:,:,k));

        %Update dFdu
        fhand=@(t,u)uotr_ctrl(t,x(k,:)',[u;modelo_u0(2:end)],1,[],modelo_par);
        fval=fhand([],uf(1));
        if imag(fval)>=1e-6
            error('Solution of states is imaginary')
        end                            
        [dFdu(:,:,k),~]=numjac(fhand,[],uf(1),real(fval),zeros(Nu,1));
        dFdu(:,:,k)=real(dFdu(:,:,k));
    else
        %No updates
        dFdx(:,:,k)=dFdx(:,:,k-1);
        dFdu(:,:,k)=dFdu(:,:,k-1);
    end
end
%}
end