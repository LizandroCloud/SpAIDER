function [r_ix2a yi] = g_ix2a(ix2,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is,in_SPAIDER)
yi = [];
C_x=1;
u_p=ones(1,1000);
figure(ix2);
tspan = linspace(t0,tf_p,1000);
ts_x(1:ns(ic)) = ts(ic,2:ns(ic)+1);
C_x(1:length(C)) = C(1:length(C));
   for i=1:length(tspan)
        u_p(i) = u_reg1(tspan(i),ts(ic,:),wopt(ic,:),is,ns(ic));
%         yi = interp1(ts_x,[wopt],in_SPAIDER.problem.interpolation, 'pp');
%         u(i)= ppval(yi,tspan(i));
        yi(i)=u_p(i);
   end  
 
 subplot(2,1,1);
 
 if strcmp(in_SPAIDER.problem.interpolation,'stage' )
    stairs(tspan*talf,u_p,'k'); hold on;
%     maxlev = wmaxlev(length(u_p),'db1');
%     [C,L] = wavedec(u_p,maxlev,'db1'); 
%     [w01,Cth,L,s,threshold] = spaider_wden(u_p,in_SPAIDER.wavelets.tptr,in_SPAIDER.wavelets.sorh,in_SPAIDER.wavelets.scal,maxlev,in_SPAIDER.wavelets.wname);
%     stairs(tspan*talf,w01,'--'); hold on;
 else
    for i=1:length(tspan)
        yi = interp1(ts_x,[wopt],'nearest', 'pp');
        u(i)= ppval(yi,tspan(i));
    end 
    plot(tspan*talf,u,'--'); hold on;
%     maxlev = wmaxlev(length(u_p),'db1');
%     [C,L] = wavedec(u_p,maxlev,'db1'); 
%     [w01,Cth,L,s,threshold] = spaider_wden(u,in_SPAIDER.wavelets.tptr,in_SPAIDER.wavelets.sorh,in_SPAIDER.wavelets.scal,maxlev,in_SPAIDER.wavelets.wname);
%     plot(tspan*talf,w01,'--'); 
 end
 
 hold on;
 figure(ix2);
  subplot(2,1,1); 
 stem(ts_x*talf,wopt(ic,1:ns(ic)),'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
 hold on;
 dU = ( max(wopt(ic,:))- min(wopt(ic,:)) )*0.2;
 if abs(min((wopt(ic,:))))< abs(max(wopt(ic,:)))
 axis([0 tf_p*talf min(wopt(ic,:))*0.8-dU  max(wopt(ic,:))*1.2+dU ])
 else
     if min((wopt(ic,:))) > 0
        axis([0 tf_p*talf  min(wopt(ic,:))*0.8-dU   max(wopt(ic,:))*1.2+dU ] )
     else
%           axis([0 tf_p*talf  min(wopt(ic,:))*0.8-dU   max(wopt(ic,:))*1.2+dU ] )
%         axis([0 tf_p*talf  max(wopt(ic,:))*1.2+dU   min(wopt(ic,:))*0.8-dU ] )   
     end
 end
 
 if nUt==1
 title(['Dynamic control profile for: ',int2str(ns_a(ic,:)),' stages']); 
 else
     title(['Dynamic control profile: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
 end
 ylabel('Control variable')
 xlabel('Time') 
 legend('stages','discrete points');
r_ix2a=1;
end