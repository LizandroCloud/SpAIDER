function[is] =  plotp(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,opC,adap,frozen,solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,wx,tx,n_est,Fobj,time,time_ac,wfact,tunn,in_SPAIDER)
  % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www. ????
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 %                              Description
 % This code is applyed to solve dynamic optimization problems with
 % wavelets basis as adaptive startegy for mesh discretization. In 
 % this strategy wavelets are used to identify prospective locations
 % at control vector parameterization where de optimization algorithm
 % will refine the discretization. To do so, details coefficents of
 % wavelets are thresholded according to a specific metric, and the
 % thesholded coefficients indicates these locations. For more information,
 % please read the article:
 % http://www.sciencedirect.com/science/article/pii/B9780444595201501032
 % ************************************************************************
 
 % changes:
 % 07-27-2012: organization and generalisation for nU control variables...
 % 08-14-2012: incorporation of tx and w...
 % 09-20-2012: improving meshref function: consideration of bounds effects
 % in refinement...
 %             improving the wavegraf function: incorporation of
 %             approximations and details pictures...
  % 11-25-2012: incorporporation of dyn.m and scipt.m for automatic derivation. 

 
%  %figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
%  adapt=1;
 w01=[];
 w02=[];
 li=0;
 ci=0;
 fref=1e5;
 w0ns=[];
 j =0;
 ns=ones(nU,1);
 wL=ones(nU,1);
 wU=ones(nU,1);
 w0=ones(nU,1);

 wLL=ones(1,nU*1);
 wUU=ones(1,nU*1);
 wopti=ones(nU,1);
 wxi=ones(nU,1);
 w=ones(nU,1);
 

 ts=ones(1,1,1);
 tss=ones(1,nU*1);

% w(1:2^isF,1:isF)=0;
% tx(nU,2^isF,10)=0;
% tol_end = 100;
kit = 1;
iter=[1 2 3 4 5 6 7 8 9 10];
% Fobj=zeros(1,10);
% time=zeros(1,10);
   
%  time_ac=zeros(1,10);
 time=zeros(1,10);
 plus=0;
%  for i=2:isF
%   time_ac(i)=time(i)+plus;
%   plus=time_ac(i);
%  end

MSE=zeros(1,10);
SNR=zeros(1,10);
toc_ant=0; 

conv =10;




 if adap == 1
     file = 'wavelets-visu';
 else
     file = 'heuristic';
 end
 % tunning of wavelet procedure...
 if opC==1  % case wavelets procedure
     if adap==1  % case adaptation on
         valK = 0.8;  % tunning parameter (valk: increase the correction of wavelet threshold)
     else % case adaptation off
         valK = 0.8;

     end
 else  % case no wavelets 
     valK = 0.0;
 end

%  diary(file)
%************************************************************************
 
thrs = []; % matrix of calculated thresholds 








 




for is = is0:isF-1 % number of stages is ns^is0 



    if is==is0   % Initial (first guess) configuration...
        for j = 1:nU; % for each control, variable
        ns(j) = (2^is); % initial number of stages
        wL(j,1:ns(j)) = u_lb(j)*ones(ns(j),1);  % lower boundary of u
        wU(j,1:ns(j)) = u_ub(j)*ones(ns(j),1);  % upper boundary of u
        w0(j,1:ns(j)) = u0(j)*ones(ns(j),1);  % initial estimative da u...
        ts(j,1:ns(j)+1) = (t0:(tf-t0)/ns(j):tf); % initial time discretization...
        w_frozen = w0;
        check = [2 3 4 5 6 7 8 9];
        ts_frozen = ts;
        ns_a = ns; 
        ts_a = ts;
        end
    else
        ts=[];
        wopt=[];
        tss=[];
    end
%  iu=is+2; % graphical counter....

for i = 1:nU
wopt(i,1:n_est(i,is-1)) = wx(i,1:n_est(i,is-1),is);
ts(i,1:n_est(i,is-1)+1) = tx(i,1:n_est(i,is-1)+1,is);
ns_a(i) = n_est(i,is-1);
ns(i)=ns_a(i);
end

disp(sprintf('Iteracao: %d', is));
disp(sprintf('Numero de estagios: %d', ns));
disp(sprintf('Método: %s',file)) 
 
 % calling optimization algorithm...
 tic % time counter
 k=0;
 tss(1) = ts(1,1);
 ki=1;
 w00=ones(1,nU*1);
 
 for i=1:nU
     
     for j=k+1:ns_a(i)+k
         tss(j)=ts(i,ki);
         ki=ki+1;
     end

     tss(j+1)=ts(i,ki);
     k = j+1;
     ki=1;
 end
 k=0; ki=1;
 
 for i=1:nU
     for j=k+1:ns(i)+k
         w00(j)=w0(i,ki);
         ki=ki+1;
     end
     k = j;
     ki=1;
 end
 
   k=0; ki=1;
    for i=1:nU
     for j=k+1:ns(i)+k
         wLL(j)=wL(i,ki);
         wUU(j)=wU(i,ki)';
         ki=ki+1;
     end
     k = j;
     ki=1;
    end
 ts_a = ts_frozen;
 nss = sum(ns);
 











 for i = 1:nU
wopt(i,1:n_est(i,is-1)) = wx(i,1:n_est(i,is-1),is);
ts(i,1:n_est(i,is-1)+1) = tx(i,1:n_est(i,is-1)+1,is);
































































end




% [wopt, erropt, iout,output,lambda,grad,hessian] = fmincon( @(ws)objfunc(x0,nss,tss,ws,ns,length(x0),nU,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE), w00, [], [], [],...
% [], wLL, wUU, @(ws)constr(x0,nss,tss,ws,ns,length(x0),nU,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE), optNLP);    












%############ reconfigurating wopt for frozen option 
    if frozen==1
    kq=1;
    for i=1:nU
            j=1:ns(i);
            wopti(i,j)=wopt(kq);
            kq=kq+1;
    end
    
    ko=0;
    if is == is0
     for i=1:nU
        for j=1:ns_a(i)
%         for k = 1:length(check(i,:))
%             if j == check(i,j)
%                   wx(i,j)=ws(ko);
%                   ko=ko+1;
%             else
                  wxi(i,j)=wopt(ko+1);
                  ko=ko+1;
%             end
        end
        
    end
else
    for i=1:nU
            for j=1:ns_a(i)
                    for k = 1:length(check(i,:))
                        if j == check(i,k)
                            wxi(i,j)=wopti(i,k);
                            break;
                        else
                            wxi(i,j)=w_frozen(i,j);
                        end
                    end
             end
    end
    end
    
%     wopt=ones(1,nU*1);
    ts=ts_a;
    wopt=wxi;
    ns=ns_a;
else
kq=1;
for i=1:nU
    for j=1:ns(i)
        wxi(i,j)=wopt(kq);
        kq=kq+1;
    end
end
wopt=[];
wopt=wxi;
ts_a=ts;
end;
for i = 1:nU
wopt(i,1:n_est(i,is-1)) = wx(i,1:n_est(i,is-1),is);
ts(i,1:n_est(i,is-1)+1) = tx(i,1:n_est(i,is-1)+1,is);
end

%#########################################################
    ic = 1:nU;
%     w(ic,1:length(wopt(ic,:)),is)=wopt(ic,:);  % each line of w is a control variable...
%     tx(ic,1:length(ts(ic,:)),is)=ts(ic,:);

 
%  toc % time counter
 %wopt = wopt'; % actual optimal control profile...
 % variable storage for future ploting... 
%     koz=1;
%     for i=1:length(ns)
%     for j=1:ns(i)
%     wopt_n(i,j)=wopt(koz);
%     koz=koz+1;
%     end    
%     end
%     wopt = [];
%     wopt = wopt_n;
 for ic = 1:nU
 
 % ploting results...
 %[ureg1 , ureg0 ] = plotfunc( x0, ns, ts, wopt, optODE,iu ); 

if opC=='on' % Wavelets analisys

if (is>is0) && (is<isF)    % if the iteration is superior than first iteration 
wopt(ic,1:length(ns(ic))) = wopt(ic,1:length(ns(ic)))';
%#########################################################################
% Rigorous Wavelets Analysis of Control Variable ....

% Solving in wavelets domain and ...
maxlev = wmaxlev(length(wopt(ic,:)),wname); % maximum level of wavelet resolution 
[C,L] = wavedec(wopt(ic,:),maxlev,wname);  % C are wavelets coefficients
detais0(ic,1:length(C),is-1)=C;  % details before thresholding

if plot ==1
ix = (is-is0)*10+ic-1;  % index for pictures generation
ix2 = ix*10;
ix3 = ix*100;
ix4 = ix*1000;
ts_x=[];
C_x=[];
% #########################################################################
% plotting wavelets coefficients...
% see help wavedec
        for i=1:length(C)
            ci(1,i)=C(i);
        end
        for i=1:length(L)
             li(1,i)=L(i);
        end
cfd = zeros(maxlev,length(wopt(ic,:))); 
for k = 1:maxlev % extracting wavelets details for plotting
     d = detcoef(ci,li,k); 
     d = d(ones(1,2^k),:);   
    cfd(k,:) = wkeep(d(:)',length(wopt(ic,:))); 
end 

cfd = cfd(:); 
I = find(abs(cfd)<sqrt(eps)); 
cfd(I)=zeros(size(I)); 
cfd = reshape(cfd,maxlev,length(wopt(ic,:)));

% plotting wavelets coefficients mesh

% scrsz = get(0,'ScreenSize');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(ix)
subplot(2,1,1), colormap(1-gray(128)); 
imagesc(wcodemat(flipud(cfd(1:maxlev,1:ns(ic)))));
% grid on;
% img = imagesc((wcodemat(cfd,128,'r'))); 
% set(get(img,'parent'),'YtickLabel',[]); 

 if nU==1
 title(['Original wavelets coefficients for: ',int2str(ns_a(ic,:)),' stages']); 
 else
     title(['Original wavelets coefficients - control variable: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
 end
ylabel('Resolution levels')
xlabel('Mesh length - Control Stages')
% colorbar
% hcb = colorbar('YTickLabel',...
% {'Low details','','','',...
% '','High details','',''});
% set(hcb,'YTickMode','manual')
% hold on;
% plotting wavelets coefficients curves

tspan = linspace(t0,tf_p,100);
ts_x(1:ns(ic)) = ts(ic,2:ns(ic)+1);
C_x(1:length(C)) = C(1:length(C));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure(ix3);
 subplot(2,1,1); 
 
    for i=1:length(tspan)
        u_p(i) = u_reg1(tspan(i),ts(ic,:),wopt(ic,:),is,ns(ic));
    end
 
 stairs(tspan,u_p,'k');
 hold on;   
 figure(ix3);
  
 subplot(2,1,1); 
 
%  hold on;
 stem(ts_x,wopt(ic,1:ns(ic)),'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
 
 hold on;
 
%  set(h,'Interpreter','none')
 axis([0 tf_p min(wopt(ic,:))*0.8  max(wopt(ic,:))*1.2 ])

  if nU==1
title(['Dynamic control profile for ',int2str(ns_a(ic,:)),' stages']); 
 ylabel('Control variable')
 xlabel('Time')
 else
     title(['Dynamic control profile: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
     ay = strcat('Control variable-',int2str(ic));
     ylabel(ay)
     xlabel('Time')
  end


 legend('stages','discrete points',1);
%  hold on;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(ix2);
    for i=1:length(tspan)
        u_p(i) = u_reg1(tspan(i),ts(ic,:),wopt(ic,:),is,ns(ic));
    end
    
 subplot(2,1,1); 
 stairs(tspan,u_p,'k');
 hold on;
 figure(ix2);
  subplot(2,1,1); 
 stem(ts_x,wopt(ic,1:ns(ic)),'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
 hold on;
 

%  set(h,'Interpreter','none')
 axis([0 tf_p min(wopt(ic,:))*0.8  max(wopt(ic,:))*1.2 ])
 
 if nU==1
 title(['Dynamic control profile for: ',int2str(ns_a(ic,:)),' stages']); 
  ylabel('Control variable')
 xlabel('Time') 
 else
     title(['Dynamic control profile: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
 
     ay = strcat('Control variable-',int2str(ic));
     ylabel(ay)
      xlabel('Time')
 end

 legend('stages','discrete points',1);
 %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tx0=[];
ko=0;
figure(ix2);
subplot(2,1,2);
tss = [1:length(C_x)];
axis([0 length(C) 1.2*min(C)  max(C)*1.2 ]);

stairs(tss,zeros(1,length(C_x)),'w');
% legend('Ap')
hold on;

for i=1:(length(L)-1)
    
     figure(ix2);
     subplot(2,1,2); 
      axis([0 length(C)+1 1.2*min(C)  max(C)*1.2 ])

    if i==1 
%         legend('Ap','')
        barcolor = 'k'; 
        elseif i==2 barcolor = 'b'; 
        legend('Ap','Dt.1')
        elseif i==3 barcolor = 'g';
        legend('Ap','Dt.1','Dt.2')
        elseif i==4 barcolor = 'r'; 
        legend('Ap','Dt.1','Dt.2','Dt.3')
        elseif i==5 barcolor = 'y';
        legend('Ap','Dt.1','Dt.2','Dt.3','Dt.4')
        elseif i==6 barcolor = 'm';   
        legend('Ap','Dt.1','Dt.2','Dt.3','Dt.4','Dt.5')
        elseif i==7 barcolor = 'c';   
    end
    
    ts_x1=(ko+1:ko+L(i));
    C_x1 = ( C_x(ko+1:ko+L(i)) );
    if length(ts_x1) == 1
        ts_x1 = [ts_x1 ( ko+ L(i)+1)];
        C_x1 = [(C_x1) (C_x( ko+ L(i)+1))];
    end
    
    bar(ts_x1,C_x1,0.8,'stacked',barcolor); 
    axis([0 length(C)+1 1.2*min(C)  max(C)*1.2 ])
    ko = ko+L(i);
    hold on;


end
grid on;
i=length(L)-1;
if i==1 
    %    legend('Ap')
       
        elseif i==2 
        legend('Legend','Ap','Dt.1')
        elseif i==3 
        legend('Legend','Ap','Dt.1','Dt.2')
        elseif i==4 
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3')
        elseif i==5
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3','Dt.4')
        elseif i==6 
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3','Dt.4','Dt.5')
        elseif i==7 
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3','Dt.4','Dt.5','Dt.6')
    end

% axis([0 tf min(C_x) max(C_x) ])

if nU==1
 title('Wavelets aproximation and details coefficients'); 
 else
     title(['Wavelets aproximation and details coefficients  - control variable: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
 end
ylabel('Wavelets coefficients')
xlabel('Points')

% subplot(3,1,3);
% bar3(ts_x,C_x); 
% axis([0 tf min(C_x)*1.2 max(C_x)*1.2 ])
% title('3D Wavelets details coefficients'); 
% ylabel('Wavelets coefficients')
% xlabel('Time')

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
 tx0=[];
ko=0;
figure(ix3);
subplot(2,1,2);
tss = [1:length(C_x)];
axis([0 length(C) 1.2*min(C)  max(C)*1.2 ]);

stairs(tss,zeros(1,length(C_x)),'w');
% legend('Ap')
hold on;
for i=1:(length(L)-1)
    
     figure(ix3);
     subplot(2,1,2); 
      axis([0 length(C)+1 1.2*min(C)  max(C)*1.2 ])
   
       if i==1 
%         legend('Ap','')
        barcolor = 'k'; 
        elseif i==2 barcolor = 'b'; 
%         legend('Ap','Dt.1')
        elseif i==3 barcolor = 'g';
%         legend('Ap','Dt.1','Dt.2')
        elseif i==4 barcolor = 'r'; 
%         legend('Ap','Dt.1','Dt.2','Dt.3')
        elseif i==5 barcolor = 'y';
%         legend('Ap','Dt.1','Dt.2','Dt.3','Dt.4')
        elseif i==6 barcolor = 'm';   
%         legend('Ap','Dt.1','Dt.2','Dt.3','Dt.4','Dt.5')
        elseif i==7 barcolor = 'c'; 
    end
     
    ts_x1=(ko+1:ko+L(i));
    C_x1 = C_x(ko+1:ko+L(i));
    if length(ts_x1) == 1
        ts_x1 = [ts_x1 ( ko+ L(i)+1)];
        C_x1 = [C_x1 C_x( ko+ L(i)+1)];
    end
    
    bar(ts_x1,C_x1,0.8,'stacked',barcolor); 
    axis([0 length(C)+1 1.5*-max(C)  1.5*max(C)])
    ko = ko+L(i);
    hold on;

end
grid on;
i=length(L)-1;
    if i==1 
    %    legend('Ap')
       
        elseif i==2 
        legend('Legend','Ap','Dt.1')
        elseif i==3 
        legend('Legend','Ap','Dt.1','Dt.2')
        elseif i==4 
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3')
        elseif i==5
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3','Dt.4')
        elseif i==6 
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3','Dt.4','Dt.5')
    end

% axis([0 tf min(C_x) max(C_x) ])
 
if nU==1
 title('Threshold values');
 else
      title(['Threshold values  - control variable: ',int2str(ic),' for ',int2str(ns_a(ic,:)),' stages'])
 end
ylabel('Wavelets coefficients')
xlabel('Points')

% subplot(3,1,3);
% bar3(ts_x,C_x); 
% axis([0 tf min(C_x)*1.2 max(C_x)*1.2 ])
% title('3D Wavelets details coefficients'); 
% ylabel('Wavelets coefficients')
% xlabel('Time')

% #########################################################################
end
% #########################################################################

j=1;
 
if adap==0 % adap =0 (fixed)
 Nc = norm(C(L(2)+1:length(C)))/length(C(L(2)+1:length(C)));  % coefficient norm...
%    Nc = mean(C(L(2)+1:length(C)));
%       Nc = 1e-6;
%     [Cth Nc] = mingcvhard(C(2:end));
  thr = ones(1,maxlev)*Nc ; % threshold value...
    thrs(ic,1:length(thr),is-is0)=thr;
        for i=L(1)+1:length(C)
            if abs(C(i))<Nc;  %fixed threshold criteria
                C(i)=0.0; % forcing these details coefficcients to be zero
            end  
        end
w01(ic,:) = waverec(C,L,wname); % inverse wavelet transformation...

% [Cth thr] = mingcvsoft(C);



else %adp=1     %(adaptive)
    
     [Cth, thrs, detais1, difdd, w01, s, threshold] = wav_thresh(C,L,maxlev,wname,detais0,wopt,...
                tptr,sorh,adap,scal,ic,is,is0,in_SPAIDER.wavelets.norm_threshold,in_SPAIDER,ns); % Wavelets thresholding
     
    
    
%      [w01(ic,:),C,L,s,threshold] = spaider_wden(wopt(ic,:),tptr,sorh,scal,maxlev,wname); % denoising
%        threshold=threshold(end:-1:1)
%      thr = thselect(wopt(ic,:),tptr)*ones(1,maxlev).*s2;  % threshold value...
        thrs(ic,1:length(thrs),is-is0)=thrs;  % threshold vector
end
%     wth(1:length( w0(ic,1:length( ns (ic) ) )),is-1)=w01(ic,1:length( ns (ic) ) );  % thresholded control profile


if plot==1
   
    figure(ix4);
    for i=1:length(tspan)
        u_pp(i) = u_reg1(tspan(i),ts(ic,:),w01(ic,:),is,ns(ic));
    end
 fit(is,:,ic) =   u_pp;
 subplot(2,1,1); 
 stairs(tspan,u_pp,'k');
 hold on;
 figure(ix4);
 subplot(2,1,1); 
 stem(ts_x,w01(ic,1:ns(ic)),'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
 legend('stages','discrete points',1);
%  grid on;
%  set(h,'Interpreter','none')
 axis([0 tf_p min(wopt(ic,1:ns(ic)))*0.8  max(wopt(ic,:))*1.2 ])
 
 if nU==1
  title(['Thresholded dynamic control profile for: ',int2str(ns_a(ic,:)),' stages']) 
   ylabel('Control variable')
 xlabel('Time') 
 else
     title(['Thresholded dynamic control profile: ',int2str(ic), ' for ',int2str(ns_a(ic,:)),' stages'])
     ay = strcat('Control variable-',int2str(ic));
     ylabel(ay)
      xlabel('Time')
 end

 
end




for i=L(2):length(C)
  ci(1,i)=C(i);
end
 for i=1:length(L)
  li(1,i)=L(i);
 end
detais1(ic,1:length(C),is-1)=C;

difd(ic,1:length(C),is-1)= detais1(ic,1:length(C),is-1)-detais0(ic,1:length(C),is-1);




% cfd = zeros(maxlev,length(wopt(ic,: ))); 
            [Pgrad] = adapt(w01,C,L,maxlev,ic,tunn,ts,adap); % prospective points...
%             Pgrad = Pgr(:,is,ic);
                        Pgrad=Pgrad(Pgrad~=0)
            PgradX(ic)=norm(Pgrad);
            if ( isempty(PgradX(ic))) 
                 elseif PgradX(ic)==0
                 else
              [Pgrad] = npoints(Pgrad,wopt(ic,:));
            end
            
%             if output.iterations < 10 
%                 Pgrad = downsample(Pgrad,2)
%             end
%           stop barrier
%             PgradX(ic)=norm(Pgrad);
            if ( isempty(PgradX(ic))) 
                stop_flag(ic)=0;
            elseif PgradX(ic)==0
                stop_flag(ic)=0;
            else
                stop_flag(ic)=1;
            end 
            if ic==nU
                stop_flag_norm=norm(stop_flag);
                if ( isempty( stop_flag_norm )) 
                     avc=1
                    return
                elseif   stop_flag_norm ==0
                    avc=1
                    return
                end
            end
            
w02(ic,: ) = wopt(ic,: ); % backing to original optimal profile...
if plot==1
% #########################################################################
% plotting wavelets coefficients...
cfd = zeros(maxlev,length(wopt(ic,:))); 
for k = 1:maxlev % extracting wavelets details for plotting
     d = detcoef(ci,li,k); 
     d = d(ones(1,2^k),:);   
    cfd(k,:) = wkeep(d(:)',length(wopt(ic,:))); 
end 

cfd = cfd(:); 
I = find(abs(cfd)<sqrt(eps)); 
cfd(I)=zeros(size(I)); 
cfd = reshape(cfd,maxlev,length(wopt(ic,:)));

% plotting wavelets coefficients
%%
figure(ix)
subplot(2,1,2), colormap(1-gray(128)); 
imagesc(wcodemat(flipud(cfd(1:maxlev,1:ns(ic)))));
%img = image(flipud(wcodemat(cfd,128,'r'))); 
%set(get(img,'parent'),'YtickLabel',[]); 

if nU==1
 title('Thresholded wavelets coefficients') 
 else
      title(['Thresholded wavelets coefficients  - control variable: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
 end
ylabel('Resolution level')
xlabel('Mesh length - Control Stages')
% colorbar
hold on;

figure(ix4)
subplot(2,1,2) 
 hist(difd(ic,1:length(C),is-1));   
 
 if nU==1
 title('Distribution of details smaller than threshold'); 
 else
    title(['Distribution of details smaller than threshold  - control variable: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
 end
 ylabel('Values of details')
 xlabel('Frequence')    



end




%#########################################################################

  % preparing for mesh refinement...
  
else  % se for menor que is0
     w02(ic,: )=wopt(ic,: ) ;
%     w0ns(ic,: )=wopt(ic,: );
     
     Pgrad = [ 1 2 3 4 5 6 7 8 ];   % Initial refinement
             
end

else % no wavelets
     if is>is0  % for first iteration...
        for i = 1:length(wopt(ic,: ))
            wopt(ic,i) = wopt(ic,i) + rand(1)*1e-4;  % forcing swmall variability...
            w02(ic,:)=wopt(ic,: );
        end;
     else
          w02(ic,1:ns(ic))=wopt(ic,1:ns(ic));
         w0ns(ic,1:ns(ic) )=wopt(ic,1:ns(ic));
     end

end

jj=1;

jj=1;

% grid refinement  (meshref function)...
%##########################################################################
[ts_ref w_ref ns_ref point] = meshref2(Pgrad,ns,ic,ts_a,w02,t0,tf,wfact);

ts(ic,1:ns_ref+1)=ts_ref;
w0ns(ic,1:ns_ref)=w_ref;
ns(ic)=ns_ref;  
 




 tam = 0;
 for i=1:length(ts(ic,:))  % constrtucting the new mesh (multivariable control)
     diff = tf-ts(ic,i);
     if (diff) <= 0.001
         tam = i;
         ns(ic) = tam-1;
         break;
     end 
 ns(ic) = 2;
 
 end
 
if (is>is0) && (is<isF) 
    %%
if plot ==1
figure(ix3);
subplot(2,1,1)
ts_x=[];

tspan = linspace(t0,tf_p,100);
ts_x(1:ns(ic)) = ts(ic,2:ns(ic)+1);


    for i=1:length(tspan)
        u_p(i) = u_reg1(tspan(i),ts(ic,:),w0ns(ic,:),is,ns(ic));
    end
    
 for i=1:ns(ic,:)
     if point(i)==1
        w0ns3(ic,i)=w0ns(ic,i);
     else
        w0ns3(ic,i)=0;
     end
 end
 stem(ts_x,w0ns3(ic,1:ns(ic)),'LineStyle',':','Marker','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
%  legend('Stages','Discret points','Refinement region','',1);
 hold on;
 stem(ts_x,w0ns(ic,1:ns(ic)),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w', 'MarkerSize',4); 
 legend('Stages','Discret points','Refinement region','New points',1);
%  set(h,'Interpreter','none')
 axis([0 tf_p min(wopt(ic,:))*0.8  max(wopt(ic,:))*1.2 ])
 if nU==1
 title(['Refined dynamic control profile for: ',int2str(ns_a(ic,:)),' stages']);
 ylabel('Refined control variable')
 xlabel('Time')
 else
      title(['Refined dynamic control profile: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
      ay = strcat('Refined control variable-',int2str(ic));
      ylabel(ay)
       xlabel('Time')
 end

hold on
%%
figure(ix3);
 subplot(2,1,1); 
 
    for i=1:length(tspan)
        u_p(i) = u_reg1(tspan(i),ts(ic,:),wopt(ic,:),is,ns_a(ic));
    end
 
 u_opt(is,:,ic) = u_p; 
%  stairs(tspan,u_p,'k');
 hold on;   
 figure(ix3);
 subplot(2,1,1); 
%  subplot(2,1,1); 
%  stairs(tspan,u_p,'k');
% %  hold on;
 stem(ts_a(ic,2:ns_a(ic)+1),wopt(ic,1:ns_a(ic)),'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
 
 
%  legend('stages','discrete points',1);
%  set(h,'Interpreter','none')
 axis([0 tf_p min(wopt(ic,:))*0.8  max(wopt(ic,:))*1.2 ])
 if nU==1
 title(['Refined dynamic control profile for: ',int2str(ns_a(ic,:)),' stages']);
  ylabel('Control variable')
 xlabel('Time')
 else
     title(['Refined dynamic control profile ',int2str(ic), ' for: ',int2str(ns_a(ic,:)),' stages'])
     ay = strcat('Refined Control variable-',int2str(ic));
     ylabel(ay)
     xlabel('Time')
 end

%  hold on;

%%
subplot(2,1,2)
ko=1;
barcolor = 'w';
for i=2:(length(L)-1)
    
    ts_x1=(ko+1:ko+L(i));
    %C_x1 = C_x(ko+1:ko+L(i));
%     if length(ts_x1) == 1
%         ts_x1 = [ts_x1 ( ko+ L(i)+1)];
%         C_x1 = [C_x1 C_x( ko+ L(i)+1)];
%     end
    %bar(ts_x1,C_x1,0.8,'stacked',barcolor); 
     if (i>1)
%         
%          if i==1 
% % %         legend('Ap','')
%          barcolor = 'k'; 
% %         elseif i==2 barcolor = 'b'; 
% %          legend('Thesh App.','Thesh Lev 1','Thesh Lev 2')
%          elseif i==3 barcolor = 'g';
% %          legend('Legend','Thesh Lev 1','Thesh Lev 2','Thesh Lev 3')
%          elseif i==4 barcolor = 'r'; 
% %         legend('Legend','Thesh Lev 0','Thesh Lev 1','Thesh Lev 2','Thesh Lev 3')
%          elseif i==5 barcolor = 'y';
% %         legend('Legend','Thesh Lev 0','Thesh Lev 1','Thesh Lev 2','Thesh Lev 3','Thesh Lev 4')
%          elseif i==6 barcolor = 'm';   
% %         egend('Legend','Thesh Lev 0','Thesh Lev 1','Thesh Lev 2','Thesh Lev 3','Thesh Lev 4','Thesh Lev 5')
%      end
        
        
    theshold=ones(1,L(i))*thrs(ic,i-1,is-is0);
    bar(ts_x1',theshold,0.5,barcolor);
%     hold on;
    bar(ts_x1',-(theshold),0.5,barcolor);
   
    ko = ko+L(i);
    
%     if i==1 
%     %    legend('Ap')
%        
%         elseif i==2 
%         legend('Legend','Thesh Lev 1')
%         elseif i==3 
%         legend('Legend','Thesh Lev 1','Thesh Lev 2')
%         elseif i==4 
%         legend('Legend','Thesh Lev 1','Thesh Lev 2','Thesh Lev 3')
%         elseif i==5
%         legend('Legend','Thesh Lev 1','Thesh Lev 2','Thesh Lev 3','Thesh Lev 4')
%         elseif i==6 
%         legend('Legend','Thesh Lev 1','Thesh Lev 2','Thesh Lev 3','Thesh Lev 4','Thesh Lev 5')
%     end
    
    
    end
    
    
    hold on;
end
  axis([0 length(C)+1 -max(C)*1.05  max(C)*1.05 ])
end




end
 %#########################################################################
 
 wL(ic,1:ns(ic)) = u_lb(ic)*ones(ns(ic),1);  % lower boundary of u
 
%  wL(ic,1:ns(ic)) = 0.8*w0ns(ic,1:ns(ic));  % lower boundary of u
%  wU(ic,1:ns(ic)) = max(1,1.1*w0ns(ic,1:ns(ic)));  % lower boundary of u
  wU(ic,1:ns(ic)) = u_ub(ic)*ones(ns(ic),1);  % upper boundary of u
 %w0(ic,1:ns(ic)) = w0ns(ic,:);

 wref(ic,1:length(w0(ic,1:length( ns (ic) ) )),is-1)=w02(ic,1:length( ns (ic) ) );
 
 end

% Fobj(is) = erropt;
 if frozen==1  % for frozen the points, only additional points will be used for optimization.
     check=[];
     w_frozen = [];
     ts_frozen=[];
     w_froz=[];
     wL_froz=[];
     wU_froz=[];
     t_froz = [];
     ns_a=[];
     ns_a = ns;
 for i=1:nU
     
        t_froz2=setxor(ts_a(i,:),ts(i,:));
        t_froz=setdiff(ts(i,:),ts_a(i,:));
            
            if (length( t_froz ) >  length(ts(i,:)) )
                t_froz = ts(i,:);
            end
        
        for j=1:length(t_froz)
        check(i,j) = find(t_froz(j) == ts(i,:));
        w_froz(i,j)=w0ns(i,check(i,j));
        wL_froz(i,j)=wL(i,check(i,j));
        wU_froz(i,j)=wU(i,check(i,j));
        end
        ns(i) = length(t_froz);
        t_froz = [];
 end
  w_frozen=w0ns;
  ts_frozen=ts;
  wL = []; wU = []; w0ns = [];
  wL = wL_froz;
  wU = wU_froz;
  w0ns = w_froz;
 else
    ns_a = ns; 
    w_frozen=w0ns;
    ts_frozen=ts;
 end
 
 li=0;
 ci=0;
 C=0;
 L=0;
 w01=[];
 w02=[];
% wopt=[];
 ts_a = [];
 w0 = w0ns;
 w0ns=[];
 wLL=[];
 wUU=[];
 ts_r =[];
 
     



 
end
 

 
         
         
         if plot==1 && is==isF-1
figure(1) % plot objective function at each iteration
subplot(2,1,1)
barcolor='k';

   bar(iter,abs(Fobj(1:length(iter))),0.1,barcolor); 

axis([0 10 min(abs(Fobj(is0+1:is-1)))*0.99  max(abs(Fobj))*1.01 ])
grid on;
title('Objective function x Iterations') 
ylabel('Objective function')
xlabel('Iterations')
%%
figure(1) % plot CPU time at each iteration
subplot(2,1,2)
time = time/(1.3*max(t_ref));
time_ac = time_ac/(1.3*max(t_ref));


barcolor='k';

   bar(iter,time_ac(1:length(iter)),0.1,barcolor);  
   hold on;
   bar(iter,time,0.1,'w'); 

axis([0 10 min(abs(time_ac))*0.9  max(abs(time_ac))*1.1 ])
grid on;
title('Normalized CPU time x Iterations') 
ylabel('Normalized CPU time')
xlabel('Iterations')
legend('CPU Accumulated','CPU on iteration',1);
%%  EMS calculation
k=1;
figure(2) % plot EMS at each iteration  
    barcolor='k';
    for ic=1:nU 
        for i=is0+1:is-1
               MSE(ic,i) = (norm(fit(i,:,ic)-u_opt(i,:,ic)))^2/100;
               SNR(ic,i) = 10*log10( (norm(u_opt(i,:,ic))^2)/100/MSE(ic,i) );
        end
               subplot(2,1,1)
               bar(iter,MSE(ic,:),0.1,barcolor);  
               if nU==1
               title('MSE x Iterations of variable')     
               else
               title(['MSE x Iterations of control variable',int2str(ic)]) 
               end
               ylabel('MSE')
               xlabel('Iterations')          
               axis([0 10 min( MSE(ic,(is0+1:is-1)))*0.9  max(MSE(ic,:))*1.1 ])
               grid on;
               subplot(2,1,2)
               bar(iter,SNR(ic,:),0.1,barcolor);  
               if nU==1
               title('SNR x Iterations')    
                    else
               title(['SNR x Iterations of control variable',int2str(ic)]) 
               end
               ylabel('SNR')
               xlabel('Iterations')          
               axis([0 10 min(SNR(ic,(is0+1:is-1)))*0.9  max(SNR(ic,:))*1.1 ])
               grid on;
    end

         end

 
  
 
  
  

 
 
 
