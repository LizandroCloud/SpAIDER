
function [r0 tstat fstat] = g_constraints(plot_state,nss,tss,wopt,ns,x0,nUt,...
    w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,f,t0,tf,tf_p,...
    talf,reg,est,wLL, wUU,optSPAIDER)
 % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: 
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 %                              Description
 % Function to plot state variables...
 % input:
 % plot_state: option for plotting
 % nss: vector of number of stages
 % tss: vector of stage time
 % wopt: vector of last calculated optimal control
 % ns: matrix (in case of two or more control variables) of number of
 % stages
 % x0: initial values
 % nUt: number of control variables
 % w_frozen:
 % 
 % ************************************************************************
 
 % changes:
 % 01-21-2013 - First modification
 %*************************************************************************
% if plot_state==1  % plotting only in this case
    if is>is0 -1% only for upper iterations
        % calling integration
        [f] = function_ode( x0,nss,tss,wopt,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER);
%         g = f +  rand(length(f),2)*0.01;
%         [f] = fun(  x0,nss,tss,wopt,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est);
        
        [ c, ceq, dc, dceq ] = constr( x0,nss,tss,wopt,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER  );


        for iz=1:length(x0)
            tamt(iz)=length(f(:,iz));
        end
        for iz=1:length(x0)
        tspan_space = linspace(t0,tf_p,tamt(iz));
        ixx=iz*1e5+is;% for different is
        if plot_state==1 
        figure(ixx)
        stairs(tspan_space*talf,f(:,iz))
        axis([0 tf_p*talf min(f(:,iz))*.95  max(f(:,iz))*1.05 ])
        title(['State Variable ',int2str(iz),' for: ',int2str(is)-is0,'th iteration']); 
        ylabel('State Variable')
        xlabel('Time')
        end
        
%         hold on;
%         if plot_state==1 
%         figure(ixx)
%         stairs(tspan_space*talf,g(:,iz))
%         axis([0 tf_p*talf min(g(:,iz))*.95  max(g(:,iz))*1.05 ])
%         title(['State Variable ',int2str(iz),' for: ',int2str(is)-is0,'th iteration']); 
%         ylabel('State Variable')
%         xlabel('Time')
%         end
        end
        izz =is;
        for iz=1:length(c(:,1))
            tamt(iz)=length(c(iz,:));
        end
        
        for iz=1:length(c(:,1))
        tspan_space = linspace(t0,tf_p,tamt(iz));
        ixx=(iz+izz)*1e5+is ;% for different is
        if plot_state==1 
        figure(ixx)
        stairs(tspan_space*talf,c(iz,:))
        axis([0 tf_p*talf min(c(iz,:))*.95  max(c(iz,:))*1.05 ])
        title(['Constraint  ',int2str(iz),' for: ',int2str(is)-is0,'th iteration']); 
        ylabel('Constraint ')
        xlabel('Time')
        end
        
%         hold on;
%         if plot_state==1 
%         figure(ixx)
%         stairs(tspan_space*talf,g(:,iz))
%         axis([0 tf_p*talf min(g(:,iz))*.95  max(g(:,iz))*1.05 ])
%         title(['State Variable ',int2str(iz),' for: ',int2str(is)-is0,'th iteration']); 
%         ylabel('State Variable')
%         xlabel('Time')
%         end
        end
        
        
        
    end
% end
r0=1;  % output flag
        tstat=tspan_space*talf;
        fstat=f;
end