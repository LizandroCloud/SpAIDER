%{
%(1) Dinamica de modelo nao-linear
function [x,G] = dynmdl(t,u,ut,ufinal,modelo_x0,modelo_u0,par)
Nx=size(modelo_x0,1);
Nu=size(u,2); %Nu=size(modelo_u0,1);
kspan=size(t,1);

%Control horizon dynamics
U=zeros(kspan,1); %U=ones(kspan,1)*modelo_u0';
for k=1:kspan
	U(k,1)=interpl(t(k),[u;ufinal],ut,1e-2); %for guess 1.
    %U(k,1)=interpl(t(k),[u;ufinal],ut,1e-2)*1e3; %for guess 1. Retornando ordem de grandeza de u para uso correto no modelo.
	%U(k,1)=interpl(t(k),[u;u(end)],ut,1e-2)*1e3; %for guess 2. Retornando ordem de grandeza de u para uso correto no modelo.
end
%{
%States
x=ones(kspan,1)*modelo_x0';
for k=2:kspan
	[~,x]=ode15s(@(ts,xs)uotr_ctrl(ts,xs,U(k-1,:)',1,[],par),(t(k-1):t(k)-t(k-1):t(k))',x(k-1,:)',opt);
    x(k,:)=xs(end,:);
end
%}

%States
opt=odeset('RelTol',1e-5,'AbsTol',1e-8,'Stats','off');
[~,x]=ode15s(@(ts,xs)uotr_ctrl(ts,xs,ut,[u;ufinal],modelo_u0,1,[],par),t,modelo_x0,opt);

%{
if nargout > 1
    %Step response matrix
    %G(:,:,1)=zeros(Nx,Nu);
    H(:,:,1)=zeros(Nx,Nu);
    for k=2:kspan
        fhand=@(ts,xs)uotr_ctrl(ts,xs,ut,[u;ufinal],modelo_u0,1,[],par);
        %fhand=@(ts,xs)uotr_ctrl(ts,xs,ut,[u;ufinal]*1e3,modelo_u0,1,[],par);
        fval=fhand(t(k),x(k,:)');
        [modelo_dFdx(:,:,k),~]=numjac(fhand,t(k),x(k,:)',fval,zeros(Nx,1));

        fhand=@(ts,us)uotr_ctrl(ts,x(k,:)',ut,us,modelo_u0,1,[],par);
        fval=fhand(t(k),U(k,1));
        [modelo_dFdu(:,:,k),~]=numjac(fhand,t(k),U(k,1),fval,zeros(Nu,1));

        [T(:,:,k),D(:,:,k)]=eig(modelo_dFdx(:,:,k)); B(:,:,k)=T(:,:,k)\modelo_dFdu(:,:,k);

        for j=1:Nx
            if D(j,j,k)>=-1e-6
                H(j,:,k)=(t(k)-t(1))*B(j,:,k);
                %G(j,:,k)=(t(k)-t(k-1))*B(j,:,k);
            else
                %H(j,:,k)=D(j,j)\(1-exp(D(j,j)*(t(k)-t(1))))*B(j,:);
                H(j,:,k)=(-D(j,j,k))\(1-exp(D(j,j,k)*(t(k)-t(1))))*B(j,:,k);
                %G(j,:,k)=(-D(j,j,k))\(1-exp(D(j,j,k)*(t(k)-t(1))))*B(j,:,k) - (-D(j,j,k))\(1-exp(D(j,j,k)*(t(k-1)-t(1))))*B(j,:,k);
            end
        end
        %{
        G(:,:,k)=T(:,:,k)*G(:,:,k);
        if imag(G(:,:,k))>=1e-6
            error('xlin solution is imaginary')
        else
            G(:,:,k)=real(G(:,:,k));
        end
        %}
    end

%{
    %Convolution matrix
    G(:,:,1)=zeros(Nx,Nu);
    for k=2:kspan
        %G(:,:,k)=T(:,:,k)*(H(:,:,k-1)-H(:,:,k));
        G(:,:,k)=T(:,:,k)*(H(:,:,k)-H(:,:,k-1));   
        if imag(G(:,:,k))>=1e-6
            error('xlin solution is imaginary')
        else
            G(:,:,k)=real(G(:,:,k));
        end
    end
%}
end
%}
%{
if nargout > 1
    %Step response matrix
    G(:,:,1)=zeros(Nx,Nu);
    for k=2:kspan
        fhand=@(ts,xs)uotr_ctrl(ts,xs,ut,[u;ufinal],modelo_u0,1,[],par);
        %fhand=@(ts,xs)uotr_ctrl(ts,xs,ut,[u;ufinal]*1e3,modelo_u0,1,[],par);
        fval=fhand(t(k),x(k,:)');
        [modelo_dFdx(:,:,k),~]=numjac(fhand,t(k),x(k,:)',fval,zeros(Nx,1));

        fhand=@(ts,us)uotr_ctrl(ts,x(k,:)',ut,us,modelo_u0,1,[],par);
        fval=fhand(t(k),U(k,1));
        [modelo_dFdu(:,:,k),~]=numjac(fhand,t(k),U(k,1),fval,zeros(Nu,1));

        G(:,:,k)=expm(modelo_dFdx(:,:,k)*(t(k)-t(1)))*modelo_dFdu(:,:,k);
    end
end
%}
end
%}
%}
%{
%(2) Dinamica de modelo linear via solucao analitica impulso
function [x,G] = dynmdl(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu)
%Funciona apenas quando tstep<=1 no arquivo main

Nx=size(modelo_x0,1); kspan=size(t,1);

%Control horizon dynamics and transition matrix
U=ones(kspan,1)*modelo_u0';
for k=1:kspan
	Phi(:,:,k)=expm(modelo_dFdx*(t(k)-t(1)));
	G(:,:,k)=Phi(:,:,k)*modelo_dFdu;
	U(k,1)=interpl(t(k),[u;ufinal],ut,1e-2)*1e3; %for guess 1. Retornando ordem de grandeza de u para uso correto no modelo.
	%U(k,1)=interpl(t(k),[u;u(end)],ut,1e-2)*1e3; %for guess 2. Retornando ordem de grandeza de u para uso correto no modelo.
end

%Convolution
cnvlt=zeros(kspan,Nx);
for k=2:kspan
	cnvlt(k,:) = cnvlt(k-1,:) + (G(:,:,k)*(U(k,:)'-modelo_u0))';
end

%States
x=ones(kspan,1)*modelo_x0';
for k=2:kspan
    x(k,:) = modelo_x0' + cnvlt(k,:);
end
end
%}

%(3) Dinamica de modelo linear via solucao analitica degrau
function [x,G] = dynmdl(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu)
Nx=size(modelo_x0,1);
Nu=size(u,2); %Nu=size(modelo_u0,1);
kspan=size(t,1);

%Control horizon dynamics and step response matrix
[T,D]=eig(modelo_dFdx); B=T\modelo_dFdu;
U=ones(kspan,1); %U=ones(kspan,1)*modelo_u0';
for k=1:kspan
    for j=1:Nx
        if D(j,j)>=-1e-6
            H(j,:,k)=(t(k)-t(1))*B(j,:);
        else
            %H(j,:,k)=D(j,j)\(1-exp(D(j,j)*(t(k)-t(1))))*B(j,:);
            H(j,:,k)=(-D(j,j))\(1-exp(D(j,j)*(t(k)-t(1))))*B(j,:);
        end
    end
	U(k,1)=interpl(t(k),[u;ufinal],ut,1e-2)*1e3; %for guess 1. Retornando ordem de grandeza de u para uso correto no modelo.
	%U(k,1)=interpl(t(k),[u;u(end)],ut,1e-2)*1e3; %for guess 2. Retornando ordem de grandeza de u para uso correto no modelo.
end

%Convolution
cnvlt=zeros(kspan,Nx);
G(:,:,1)=zeros(Nx,Nu);
for k=2:kspan
    %G(:,:,k)=T*(H(:,:,k-1)-H(:,:,k));
    G(:,:,k)=T*(H(:,:,k)-H(:,:,k-1));
    if imag(G(:,:,k))>=1e-6
        error('xlin solution is imaginary')
    else
        G(:,:,k)=real(G(:,:,k));
    end
    for j=1:Nx
        if D(j,j)>=-1e-6
            cnvlt(k,j) = cnvlt(k-1,j) + (G(j,:,k)*U(k-1,:)')';
        else
            cnvlt(k,:) = cnvlt(k-1,:) + (G(:,:,k)*(U(k-1,:)'-modelo_u0(1)))';
        end
    end    
    %cnvlt(k,:) = cnvlt(k-1,:) + (G(:,:,k)*(U(k-1,:)'-modelo_u0))';
end

%States
x=ones(kspan,1)*modelo_x0';
for k=2:kspan
    x(k,:) = modelo_x0' + cnvlt(k,:);
end
end
%}