function  [r_ix4c fit] = g_ix4c(ix4,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,w01,wopt,is)
%     figure(ix4);
    tspan = linspace(t0,tf_p,1000);
%     ts_x(1:ns(ic)) = ts(ic,2:ns(ic)+1);
%     C_x(1:length(C)) = C(1:length(C));
    for i=1:length(tspan)
        u_pp(i) = u_reg1(tspan(i),ts(ic,:),wopt(ic,:),is,ns(ic));
    end
 fit(is,:,ic) =   u_pp;
%  subplot(2,1,1); 
%  stairs(tspan*talf,u_pp,'k');
%  hold on;
%  figure(ix4);
%  subplot(2,1,1); 
%  stem(ts_x*talf,w01(ic,1:ns(ic)),'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
%  legend('stages','discrete points',1);
%  grid on;
%  set(h,'Interpreter','none')
%  axis([0 tf_p*talf min(wopt(ic,1:ns(ic)))*0.8  max(wopt(ic,:))*1.2 ])
 
%  if nUt==1
%  title(['Thresholded dynamic control profile for: ',int2str(ns_a(ic,:)),' stages']) 
%  else
%      title(['Thresholded dynamic control profile: ',int2str(ic), ' for ',int2str(ns_a(ic,:)),' stages'])
%  end
%  ylabel('Control variable')
%  xlabel('Time')
 r_ix4c=1;
end