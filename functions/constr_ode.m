function [ c, ceq, dc, dceq ] = constr_ode( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER  )
wsN=ws;
 % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www. ????
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 %                              Description
 % This routine is used to transfer information of variales t calculate constraints...
 %
 % changes:
 % 07-20-2013: printing algebraic variables...
 %
 % ************************************************************************
 k=0; ki=1;
%   for i=1:n_c % number of control variables
%      for j=k+1:ne(i)+k
%          wsN(j)=ws(j)*(wUU(j)-wLL(j))+wLL(j);
%          ki=ki+1;
%      end
%      k = j;
%      ki=1;
%   end 
  
 if typ==1 % if terminal free
     wsN(n_c+1)=ws(n_c+1)*(1000-100)+100;
 end
  
if nargout == 2   % numerical jacobian
	[f, t_ult,t_out,yout,n_elem] = function_ode( x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER);
	% Modificação UOTR
    [ceq c] = cr1(f, t_ult,t_out,ws,yout,n_elem,ne); % function for calculate constraints
else              % analytical jacobian
	[f,t_ult,t_out,yout,n_elem,df] = function_ode( x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER);
    
%       % correcting sensitivy from normalization...
%      k=0; ki=1;
%      for i=1:n_c % number of control variables
%         for j=k+1:ne(i)+k
%             df(:,:,j)=df(:,:,j)*(wUU(j)-wLL(j)); 
%             ki=ki+1;
%         end
%         k = j;
%         ki=1;
%      end
    % Modificação UOTR
	[ceq c dceq dc] = cr2(f, t_ult,t_out,df,ws,yout,n_elem, ne);  %function for calculate constraints and derivatives of constraints
    
end
end