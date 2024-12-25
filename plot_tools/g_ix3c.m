function [r_ix3c u_opt] = g_ix3c(ix3,ic,ns,ns_a,t0,tf_p,ts,ts_a,talf,C,nUt,wopt,w0ns,point,is)
        figure(ix3);
        subplot(2,1,1)
%         legend('Stages','Discret points','Refinement region','New points',1);
        tspan = linspace(t0,tf_p*talf,1000);
        ts_x(1:ns(ic)) = ts(ic,2:ns(ic)+1);
            for i=1:length(tspan)
                u_p(i) = u_reg1(tspan(i),ts(ic,:)*talf,w0ns(ic,:),is,ns(ic));
            end
        for i=1:ns(ic,:)
            if point(i)==1
                w0ns3(ic,i)=w0ns(ic,i);
            else
                w0ns3(ic,i)=0;
            end
        end
        dU = ( max(wopt(ic,:))- min(wopt(ic,:)) )*0.2;
        stem(ts_x*talf,w0ns3(ic,1:ns(ic)),'LineStyle',':','Marker','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); hold on;
% 		stem(ts_x*talf,(w0ns(ic,1:ns(ic)))*1.2+dU,'LineStyle',':','Marker','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); hold on;
%         stem(ts_x*talf,(w0ns(ic,1:ns(ic)))*0.8-dU,'LineStyle',':','Marker','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
        hold on;
        stem(ts_x*talf,w0ns(ic,1:ns(ic)),'LineStyle','none','Marker','x','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
        legend('Stages','Discret points','Refinement region','New points');
        stem(ts_x*talf,(w0ns(ic,1:ns(ic)))*1.2+dU,'LineStyle','none','Marker','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); hold on;
        stem(ts_x*talf,(w0ns(ic,1:ns(ic)))*0.8-dU,'LineStyle','none','Marker','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
         if min(wopt(ic,:))<max(wopt(ic,:))
%             axis([0 tf_p*talf min(wopt(ic,:))*0.8-dU  max(wopt(ic,:))*1.2+dU ])
         else
%             axis([0 tf_p*talf  max(wopt(ic,:))*1.2+dU min(wopt(ic,:))*0.8-dU] )
         end
        if nUt==1
            title(['Refined dynamic control profile for: ',int2str(ns_a(ic,:)),' stages']);
        else
            title(['Refined dynamic control profile: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
        end
        ylabel('Refined control variable')
        xlabel('Time')
        hold on
        figure(ix3);
        subplot(2,1,1); 
        for i=1:length(tspan)
            u_p(i) = u_reg1(tspan(i),ts(ic,:)*talf,wopt(ic,:),is,ns_a(ic));
        end
         dU = ( max(w0ns(ic,:))- min(w0ns(ic,:)) )*0.2;
         u_opt(is,:,ic) = u_p; 
         stem(ts_a(ic,2:ns_a(ic)+1)*talf,wopt(ic,1:ns_a(ic)),'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); hold on;
%          stem(ts_a(ic,2:ns_a(ic)+1)*talf,(w0ns(ic,1:ns(ic)))*1.2+dU,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); hold on;
%          stem(ts_a(ic,2:ns_a(ic)+1)*talf,(w0ns(ic,1:ns(ic)))*0.8-dU,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',3); 
            if min(wopt(ic,:))<max(wopt(ic,:))
%             axis([0 tf_p*talf min(wopt(ic,:))*0.8-dU  max(wopt(ic,:))*1.2+dU ])
          	else
%             axis([0 tf_p*talf  max(wopt(ic,:))*1.2+dU min(wopt(ic,:))*0.8-dU] )
            end
            if min(wopt(ic,:))<max(wopt(ic,:))
%             axis([0 tf_p*talf min(wopt(ic,:))*0.8-dU  max(wopt(ic,:))*1.2+dU ])
            else
%             axis([0 tf_p*talf  max(wopt(ic,:))*1.2+dU min(wopt(ic,:))*0.8-dU] )
            end
         if nUt==1
            title(['Refined dynamic control profile for: ',int2str(ns_a(ic,:)),' stages']);
         else
            title(['Refined dynamic control profile ',int2str(ic), ' for: ',int2str(ns_a(ic,:)),' stages'])
         end
         ylabel('Control variable')
         xlabel('Time')
         r_ix3c=1;
end