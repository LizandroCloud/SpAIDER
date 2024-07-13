function [sc conv conv2] = spaider_conv(is,is0,ic,Fobj,f_tol,f_tol2,fref,plot,iter,time,t_ref,u_opt,conv,convv,time_ac,nUt,upt)
sc=0;
% conv = [10 10];


% based on the parameterization of control variables...
 ksum = 0;
 if is==is0+1
     conv(is-2)=abs((abs(Fobj(is))-abs(Fobj(is-1)))*abs((Fobj(is)))^-1) + 5*abs(min((( abs(Fobj(is))-abs(fref) )/(abs(fref))),0)) % first convergence criterium
     
     for i=1:length(upt(1,:))
         conv2(is-2) = 0 + ksum;
         conv2 = sqrt(conv2);
         ksum = conv2(is-2);
     end

 elseif is>is0+1
          conv(is-2)=abs((abs(Fobj(is))-abs(Fobj(is-1)))*abs((Fobj(is)))^-1) + 5*abs(min((( abs(Fobj(is))-abs(fref) )/(abs(fref))),0)); % first convergence criterium
     
     for i=1:length(upt(ic,is,:))
         for ic=1:nUt;
%          conv2 = 0;
         try
         conv2(is-2) = ( upt(ic,is,i)-upt(ic,is-1,i) )^2 + ksum;
         catch 
           conv2=ones(1,30);
         
         ksum = conv2(is-2);
%         ksum = 0;
          end
         end
     end 
     conv2(is-2) = sqrt(conv2(is-2))
     conv(is-2) = 0.001*conv2(is-2) + conv(is-2)
     
 else
     conv(is-2)=10;  % default
 end
 
 if is>is0+1
     convv(is-2)=abs((abs(conv(is-2))-abs(conv(is-3)))*abs((conv(is-3)))^-1) ; % second convergence criterium
 else
     convv(is-2)=10; % default
 end
  % -----------------------------------------------------------------------
 
 if (conv(is-2)>f_tol)  % checking first convergence
     if ((convv(is-2)>f_tol2))   % checking second convergence
        sc=0;
     elseif ((convv(is-2)<=f_tol2) )
%      elseif ((convv(is-2)<=f_tol2) && (abs(Fobj(is))>=abs(fref)))
%             if plot==1
%                 
% %                 r1 = g_function(Fobj,is0,is,time,t_ref);
% %                 r2 = g_mse(nUt,is0,is,fit,ic,u_opt);
%                 
%                 figure(1) % plot objective function at each iteration
%                 subplot(2,1,1)
%                 barcolor='k';
%                 bar(iter,abs(Fobj),0.1,barcolor); 
%                 axis([0 10 min(abs(Fobj(is0+1:is)))*0.99  max(abs(Fobj))*1.01 ])
%                 grid on;
%                 title('Objective function x Iterations') 
%                 ylabel('Objective function')
%                 xlabel('Iterations')
%                 figure(1) % plot CPU time at each iteration
%                 subplot(2,1,2)
%                 time = time/(1.3*max(t_ref));
%                 time_ac = time_ac/(1.3*max(t_ref));
%                 barcolor='k';
%                 bar(iter,time_ac,0.1,barcolor);  
%                 hold on;
%                 bar(iter,time,0.1,'w'); 
%                 axis([0 10 min(abs(time_ac(is0+1:is)))*0.9  max(abs(time_ac))*1.1 ])
%                 grid on;
%                 title('Normalized CPU time x Iterations') 
%                 ylabel('Normalized CPU time')
%                 xlabel('Iterations')
%                 legend('CPU Accumulated','CPU on iteration',1);
%                 
%                 %%  EMS calculation
%                 k=1;
%                 figure(2) % plot EMS at each iteration  
%                 barcolor='k';
%                 for ic=1:nUt 
%                     for i=is0+1:is-1
%                         MSE(ic,i) = (norm(fit(i,:,ic)-u_opt(i,:,ic)))^2/100;
%                         SNR(ic,i) = 10*log10( (norm(u_opt(i,:,ic))^2)/100/MSE(ic,i) );
%                     end
%                 subplot(2,1,1)
% %               bar(iter,MSE(ic,length(iter)),0.1,barcolor);  
%                 if nUt==1
%                     title('MSE x Iterations of variable')     
%                 else
%                     title(['MSE x Iterations of control variable',int2str(ic)]) 
%                 end
%                 ylabel('MSE')
%                 xlabel('Iterations')          
% %               axis([0 10 min( MSE(ic,(is0+1:is)))*0.9  max(MSE(ic,:))*1.1 ])
%                 grid on;
%                 subplot(2,1,2)
% %               bar(iter,SNR(ic,:),0.1,barcolor);  
%                 if nUt==1
%                    title('SNR x Iterations')    
%                 else
%                    title(['SNR x Iterations of control variable',int2str(ic)])
%                 end
%                 ylabel('SNR')
%                 xlabel('Iterations')          
% %               axis([0 10 min(SNR(ic,(is0+1:is)))*0.9  max(SNR(ic,:))*1.1 ])
%                 grid on;
%                 end
%                 grid on;
%             end
            sc=0;
     end
      elseif ( ((conv(is-2)<= f_tol) &&(convv(is-2)<=f_tol2)) )
%  elseif ( ((conv(is-2)<= f_tol) &&(convv(is-2)<=f_tol2))&& (abs(Fobj(is))>=abs(fref)) )

     
     
               %% Plotting Objective Function
if plot==1
    
%                 r1 = g_function(Fobj,is0,is,time,t_ref);
%                 r2 = g_mse(nUt,is0,is,fit,ic,u_opt);
    
figure(1) % plot objective function at each iteration
subplot(2,1,1)
barcolor='k';

    bar(iter,abs(Fobj),0.1,barcolor); 

axis([0 10 min(abs(Fobj(is0+1:is)))*0.99  max(abs(Fobj))*1.01 ])
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

   bar(iter,time_ac,0.1,barcolor);  
   hold on;
   bar(iter,time,0.1,'w'); 

% axis([0 10 min(abs(time_ac(is0+1:is)))*0.9  max(abs(time_ac))*1.1 ])
grid on;
title('Normalized CPU time x Iterations') 
ylabel('Normalized CPU time')
xlabel('Iterations')
legend('CPU Accumulated','CPU on iteration');
%%  EMS calculation
% k=1;
% figure(2) % plot EMS at each iteration  
%     barcolor='k';
%     for ic=1:nUt 
%         for i=is0+1:is-1
%                MSE(ic,i) = (norm(fit(i,:,ic)-u_opt(i,:,ic)))^2/100;
%                SNR(ic,i) = 10*log10( (norm(u_opt(i,:,ic))^2)/100/MSE(ic,i) );
%         end
%                subplot(2,1,1)
%                bar(iter,MSE(ic,:),0.1,barcolor);  
%                if nUt==1
%                title('MSE x Iterations of variable')     
%                else
%                title(['MSE x Iterations of control variable',int2str(ic)]) 
%                end
%                ylabel('MSE')
%                xlabel('Iterations')          
%                axis([0 10 min( MSE(ic,(is0+1:is)))*0.9  max(MSE(ic,:))*1.1 ])
%                grid on;
%                subplot(2,1,2)
%                bar(iter,SNR(ic,:),0.1,barcolor);  
%                if nUt==1
%                title('SNR x Iterations')    
%                     else
%                title(['SNR x Iterations of control variable',int2str(ic)]) 
%                end
%                ylabel('SNR')
%                xlabel('Iterations')          
%                axis([0 10 min(SNR(ic,(is0+1:is)))*0.9  max(SNR(ic,:))*1.1 ])
%                grid on;
%     end
% 
% 
%  
%    grid on;
 
 end
             sc=1;
 end
 
end