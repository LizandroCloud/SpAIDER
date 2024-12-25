function [r_ix3a] = g_ix3a(ix3,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is,in_SPAIDER)
C_x=1;
u_p=ones(1,1000);
figure(ix3);
tspan = linspace(t0,tf_p,1000);
ts_x(1:ns(ic)) = ts(ic,2:ns(ic)+1);
C_x(1:length(C)) = C(1:length(C));
subplot(2,1,1); 
for i=1:length(tspan)
    u_p(i) = u_reg1(tspan(i),ts(ic,:),wopt(ic,:),is,ns(ic));
%     yi = interp1(ts_x,[wopt],in_SPAIDER.problem.interpolation, 'pp');
%     u(i)= ppval(yi,tspan(i));
end

 if strcmp(in_SPAIDER.problem.interpolation,'stage' )
    stairs(tspan*talf,u_p,'k');
 else
    for i=1:length(tspan)
    yi = interp1(ts_x,[wopt],'constant', 'pp');
    u(i)= ppval(yi,tspan(i));
    end 
    plot(tspan*talf,u,'--'); 
 end
hold on;   
figure(ix3);
subplot(2,1,1); 
stem(ts_x*talf,wopt(ic,1:ns(ic)),'LineStyle','none','MarkerEdgeColor','none','MarkerFaceColor','k', 'MarkerSize',2); 
hold on;
%  set(h,'Interpreter','none')
dU = ( max(wopt(ic,:))- min(wopt(ic,:)) )*0.2;
 if min(wopt(ic,:))<= max(wopt(ic,:))
 axis([0 tf_p*talf min(wopt(ic,:))*0.8-dU  max(wopt(ic,:))*1.2+dU ])
 else
   axis([0 tf_p*talf  max(wopt(ic,:))*1.2+dU min(wopt(ic,:))*0.8-dU] )
 end
if nUt==1
    title(['Dynamic control profile for ',int2str(ns_a(ic,:)),' stages']); 
else
    title(['Dynamic control profile: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
end
ylabel('Control variable')
xlabel('Time')
legend('stages','discrete points');
r_ix3a=1;
end