    for i=1:length(tspan)
        u_p(i) = u_reg1(tspan(i),ts(ic,:),wopt(ic,:),is,ns(ic));
    end
    
 subplot(2,1,1); 
 stairs(tspan*talf,u_p,'k');
 hold on;
 figure(ix2);
  subplot(2,1,1); 
 stem(ts_x*talf,wopt(ic,1:ns(ic)),'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
 hold on;
 

%  set(h,'Interpreter','none')
 axis([0 tf_p*talf min(wopt(ic,:))*0.8  max(wopt(ic,:))*1.2 ])
 
 if nUt==1
 title(['Dynamic control profile for: ',int2str(ns_a(ic,:)),' stages']); 
 else
     title(['Dynamic control profile: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
 end
 ylabel('Control variable')
 xlabel('Time') 
 legend('stages','discrete points',1);