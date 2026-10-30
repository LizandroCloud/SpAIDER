function s = modelo_eqsens_nonlin(t_out,s0,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,dFdx,dFdu)
global modelo_dFdx modelo_dFdu

kspan=size(t_out,1);

%[T,D]=eig(modelo_dFdx); B=T\modelo_dFdu*1e3;
[T,D]=eig(dFdx); B=T\dFdu*1e3;
for k=1:kspan
    for j=1:n_s
        if real(D(j,j))>=-1e-6
            H(j,:,k)=(t_out(k)-t_out(1))*B(j,:);
            Phi(j,:,k)=zeros(1,n_s); Phi(j,j,k)=1;
        else
            H(j,:,k)=(-D(j,j))\(1-exp(D(j,j)*(t_out(k)-t_out(1))))*B(j,:);
            Phi(j,:,k)=zeros(1,n_s); Phi(j,j,k)=exp(D(j,j)*(t_out(k)-t_out(1)));
        end
    end
    St(:,:,k)=real(T*H(:,:,k));
    Et(:,:,k)=real(T*Phi(:,:,k)/T);
end

%initialize sensitivity stock vector
s=[];
for j=1:n_c
    s=[s,zeros(1,n_s*ns_a(j))];
end

%Sensitivity matrices allocation in output [dx]
j0=0;
for j=1:n_c
	for is = 1:ns_a(j)
        %id=zeros(ns_a(j),1); id(is)=1;
        idt=max(find(ts(is)>=t_out));
        if is<ns_a(j)
            idtp=max(find(ts(is+1)>=t_out));
        else
            idtp=length(t_out);
        end
        for k=idt:idtp
            %uws=interpl(t_out(k),id,ts,1e-2);
            uws=1;
            s(k,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=...
                s0(1,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0) + (St(:,j,k-idt+1)*uws)';
        end
        for k=idtp:kspan
            s(k,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=...
                (Et(:,:,k-idtp+1)*s(idtp,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)')';
        end
        if idt>1
            for k=1:idt-1
                s(k,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=...
                    s0(1,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0);
            end
        end
	end 
	%j0=size(s,2);
    j0=j0+ns_a(j)*n_s;
end

%{
j0=0;
for j=1:n_c
    for is=1:ns_a(j)
        cnvlt=zeros(kspan,n_s);
        G(:,:,1)=zeros(n_s,n_c);
        id=zeros(ns_a(j),1); idp=id;
        id(is)=1;
        if is<ns_a(j)
            idp(is+1)=1;
        end
        %idt=max(find(ts(is)>=t_out));
        %if is==ns_a(j)
        %    idtp=length(t_out);
        %else
        %    idtp=max(find(ts(is+1)>=t_out));
        %end
        for k=2:kspan
            G(:,:,k)=T*(H(:,:,k)-H(:,:,k-1));
            %G(:,:,k)=T*(H(:,:,k)-H(:,:,1));
            %G(:,:,k)=T*H(:,:,k);
            if imag(G(:,:,k))>=1e-6
                error('xlin solution is imaginary')
            end
            G(:,:,k)=real(G(:,:,k));
            
            %Incidence when ws(ns_a(j)) is active
            %uws=interpl(t_out(k-1),id,ts,1e-2);
            uws=interpl(t_out(k),id,ts,1e-2);
            uwsp=interpl(t_out(k),idp,ts,1e-2);
            %uws=1;
            
            %Convolution
            cnvlt(k,:) = cnvlt(k-1,:) + (G(:,j,k)*uws)';
            %cnvlt(k,:) = (G(:,j,k)*uws)';

            %Sensitivities
            %sp(k,:) = sp(1,:) + cnvlt(k,:);
            s(k,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=cnvlt(k,:);
            %s(k,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=(G(:,j,k)*uws)';
        end
        %Stalling sensitivities into block vector
        %s(:,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=sp;
    end
    j0=size(s,2);
end
%}
end