%{
%Dinamica de linear via solucao analitica degrau usando jacobiana variável
function s = modelo_eqsens(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,dFdx)
global modelo_u0 
global modelo_dFdu

kspan=size(t_out,1);

%Control horizon dynamics and step response matrix
[T,D]=eig(dFdx); B=T\modelo_dFdu;
U=ones(kspan,1); %U=ones(kspan,1)*modelo_u0';
for k=1:kspan
    %[T,D]=eig(dFdx(:,:,k)); B=T\modelo_dFdu;
    for j=1:n_s
        if real(D(j,j))>=-1e-6
            H(j,:,k)=(t_out(k)-t_out(1))*B(j,:);
        else
            H(j,:,k)=(-D(j,j))\(1-exp(D(j,j)*(t_out(k)-t_out(1))))*B(j,:);
        end
    end
	U(k,1)=interpl(t_out(k),ws,ts,1e-2)*1e3;
end

%Convolution
j0=0;
for j=1:n_c
    for is=1:ns_a(j)
        sp(1,:)=zeros(n_s,1)';
        cnvlt=zeros(kspan,n_s);
        G(:,:,1)=zeros(n_s,n_c);
        for k=2:kspan
            id=zeros(ns_a(j),1); id(is)=1;
            uws=interpl(t_out(k-1),id,ts,1e-2);
            G(:,:,k)=uws*T*(H(:,:,k)-H(:,:,k-1));
            if imag(G(:,:,k))>=1e-6
                error('xlin solution is imaginary')
            end
            G(:,:,k)=real(G(:,:,k));
            for i=1:n_s
                if real(D(i,i))>=-1e-6
                    cnvlt(k,i) = cnvlt(k-1,i) + (G(i,:,k)*U(k-1,:)')';
                else
                    cnvlt(k,:) = cnvlt(k-1,:) + (G(:,:,k)*(U(k-1,:)'-modelo_u0(1)))';
                end
            end
        
            %Sensitivities
            sp(k,:) = sp(1,:) + cnvlt(k,:)*tal;
        end
        
        %Stalling sensitivities into block vector
        s(:,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=sp;
    end
    j0=size(s,2);
end

end
%}

%Dinamica de modelo linear via solucao analitica degrau
function s = modelo_eqsens(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
global modelo_u0 
global modelo_dFdx modelo_dFdu

kspan=size(t_out,1);

%Control horizon dynamics and step response matrix
[T,D]=eig(modelo_dFdx); B=T\modelo_dFdu;
U=ones(kspan,1); %U=ones(kspan,1)*modelo_u0';
for k=1:kspan
    for j=1:n_s
        if real(D(j,j))>=-1e-6
            H(j,:,k)=(t_out(k)-t_out(1))*B(j,:);
        else
            H(j,:,k)=(-D(j,j))\(1-exp(D(j,j)*(t_out(k)-t_out(1))))*B(j,:);
        end
    end
	U(k,1)=interpl(t_out(k),ws,ts,1e-2)*1e3;
end

%Convolution
j0=0;
for j=1:n_c
    for is=1:ns_a(j)
        sp(1,:)=zeros(n_s,1)';
        cnvlt=zeros(kspan,n_s);
        G(:,:,1)=zeros(n_s,n_c);
        for k=2:kspan
            id=zeros(ns_a(j),1); id(is)=1;
            uws=interpl(t_out(k-1),id,ts,1e-2);
            G(:,:,k)=uws*T*(H(:,:,k)-H(:,:,k-1));
            if imag(G(:,:,k))>=1e-6
                error('xlin solution is imaginary')
            end
            G(:,:,k)=real(G(:,:,k));
            for i=1:n_s
                if real(D(i,i))>=-1e-6
                    cnvlt(k,i) = cnvlt(k-1,i) + (G(i,:,k)*U(k-1,:)')';
                else
                    cnvlt(k,:) = cnvlt(k-1,:) + (G(:,:,k)*(U(k-1,:)'-modelo_u0(1)))';
                end
            end
        
            %Sensitivities
            sp(k,:) = sp(1,:) + cnvlt(k,:)*tal;
        end
        
        %Stalling sensitivities into block vector
        s(:,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=sp;
    end
    j0=size(s,2);
end

end
%}