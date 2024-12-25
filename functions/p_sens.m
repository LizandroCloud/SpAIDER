% function s = p_sens(t_out,s0,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,dFdx,dFdu)  
% function s = p_sens(t_out,s0,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,dFdx,dFdu)  
function   s = p_sens(t_out,s0,      ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg,dFdx,dFdu)
               % p_sens(t_out,s(end,:),ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,check, tal,est,k,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref,dFdx(:,:,k),dFdu(:,:,k));
                     % % global modelo_dFdx modelo_dFdu
    kspan=length(t_out);
    % 
    %initialize sensitivity stock vector
    s(1,:)=s0;
    % 
    [T,D]=eig(dFdx); 
    B=T\dFdu; 
    dt=t_out(2)-t_out(1);
    for j=1:n_s
        if real(D(j,j))>=-1e-6
            H(j,:)=dt*B(j,:);
            Phi(j,:)=zeros(1,n_s); Phi(j,j)=1;
        else
            H(j,:)=(-D(j,j))\(1-exp(D(j,j)*dt))*B(j,:);
            Phi(j,:)=zeros(1,n_s); Phi(j,j)=exp(D(j,j)*dt);
        end
    end
    St(:,:)=real(T*H(:,:));
    Et(:,:)=real(T*Phi(:,:)/T);
    % 
    j0=0;
    for j=1:n_c
        for is = 1:ns_a(j)
            id=zeros(ns_a(j),1); id(is)=1;
            for k=2:kspan
                %Incidence when ws(ns_a(j)) is active
                uws=interpl(t_out(k-1),id,ts(j,:),1e-2);
    % 
                %Sensitivities by convolution and previous state
                s(k,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0) = ...
                    (Et*s(k-1,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)')' ...
                    + St(:,j)'*uws;
            end
        end
        j0=j0+ns_a(j)*n_s;
    end
end


% 
% %Dinamica linear via solucao analitica degrau usando jacobiana variável
% function s = p_sens(t_out,s0,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg,dFdx,dFdu)
% % global modelo_dFdx modelo_dFdu
% 
% kspan=size(t_out',1);
% 
% [T,D]=eig(dFdx); B=T\dFdu;
% % B(:,1)=B(:,1)*1e3;
% for k=1:kspan
%     for j=1:n_s
%         if real(D(j,j))>=-1e-6
%             H(j,:,k)=(t_out(k)-t_out(1))*B(j,:);
%             Phi(j,:,k)=zeros(1,n_s); Phi(j,j,k)=1;
%         else
%             H(j,:,k)=(-D(j,j))\(1-exp(D(j,j)*(t_out(k)-t_out(1))))*B(j,:);
%             Phi(j,:,k)=zeros(1,n_s); Phi(j,j,k)=exp(D(j,j)*(t_out(k)-t_out(1)));
%         end
%     end
%     St(:,:,k)=real(T*H(:,:,k));
%     Et(:,:,k)=real(T*Phi(:,:,k)/T);
% end
% 
% %initialize sensitivity stock vector
% s=[];
% for j=1:n_c
%     s=[s,zeros(1,n_s*ns_a(j))];
% end
% 
% %Sensitivity matrices allocation in output [dx]
% j0=0;
% for j=1:n_c
% 	% for is = 1:ns_a(j)
%       for is = yi:ns_a(j)
%         id=zeros(ns_a(j),1); id(is)=1;
%         idt=max(find(ts(is)>=t_out));
%         idt=1;
%         if is<ns_a(j)
%             idtp=max(find(ts(is+1)>=t_out));
%         else
%             idtp=length(t_out);
%         end
%         for k=idt:idtp
%             uws=1;
%             s(k,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=...
%                 s0(1,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0) + (St(:,j,k-idt+1)*uws)';
%             if yi>1
%                 disp(size(s0)); % Exibe as dimensões de s0
%                 disp(idt); % Exibe o valor de idt
%                 disp(1:n_s*(is-1)); % Exibe o índice que está sendo acessado
%                 s(k,1:n_s*(is-1))=...
%                     (Et(:,:,k-idt+1)*s0(idt,1:n_s*(is-1))')';
%             end
%         end
%         for k=idtp:kspan
%             s(k,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=...
%                 (Et(:,:,k-idtp+1)*s(idtp,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)')';
%         end
%         if idt>1
%             for k=1:idt-1
%                 s(k,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0)=...
%                     s0(1,n_s*(is-1)+1+j0:n_s*(is-1)+n_s+j0);
%             end
%         end
% 	end 
%     j0=j0+ns_a(j)*n_s;
% end
% end