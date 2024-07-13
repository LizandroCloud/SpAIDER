 function [ J, dJ ] = objfunc( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER)
  % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www. ????
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 %                              Description
 % This routine is used to transfer information of variales ot calculate obj fucntional...
 %
 % changes:
 % 07-20-2013: printing algebraic variables...
 %
 % ************************************************************************
 
  wsN=ws;
% correcting normalization
 k=0; ki=1;
%   for i=1:n_c % number of control variables
%      for j=k+1:ne(i)+k
%          wsN(j)=ws(j)*(wUU(j)-wLL(j))+wLL(j);
%          ki=ki+1;
%      end
%      k = j;
%      ki=1;
%   end
%  if typ==1 % if terminal free
%      wsN(n_c+1)=ws(n_c+1)*(1000-100)+100;
%  end
 if nargout == 1   % only objective functional....
 
 [f, t_ult, t_out, yout] = fun( x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER);
 
 [J]= obf1(f,t_ult,t_out,wsN,yout);

 else % objective function and constraints...
  [f,t_ult,t_out,yout,n_elem,df] = fun(  x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER);

%  cc = (1/f(end,3))*(cao*Vo - f(end,1)*f(end,3));  % calculo de cc no tempo final...
%  J = -cc*f(end,3);  
  [J dJ]= obf2(f,t_ult,t_out,df,wsN,yout); 
  
%   % correcting sensitivy from normalization...
%   k=0; ki=1;
%   for i=1:n_c % number of control variables
%      for j=k+1:ne(i)+k
%          dJ(:,:,j)=dJ(:,:,j)*(wUU(j)-wLL(j)); 
%          ki=ki+1;
%      end
%      k = j;
%      ki=1;
%   end

  
 % calculo da Fobj no tempo final...
%  dcc = ( (-1/f(3)^2)*(cao*Vo))*df(3,:) - df(1,:)*1;
% 
%  syms f3 f2 f1 cb_in cao cb_in cbo_in Vo
%  eq = [-(cao*Vo - f1*f3)];
%  [jac,j1] = jacmat( eq , [f1,f2,f3], f );
 %dJ = df(end,1,:)+ df(end,2,:);
 
 
 
 
 
 end
 end