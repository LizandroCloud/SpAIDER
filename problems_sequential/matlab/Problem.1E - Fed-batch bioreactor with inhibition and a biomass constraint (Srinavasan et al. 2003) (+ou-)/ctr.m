function [ c, ceq, dc, dceq ] = ctr( x0, ns, ts, ws, optODE )
 
 cb_in = 5;         % parametros iniciais... ( mol/l )
 Vo = 1;      
 cao = 0.72;
 cbo = 0.05;    
 cd_max = 0.15;
 cb_max = 0.025;
 
 if nargout == 2   % caso nao tenha analise de sensibilidade...
 [f] = fun( x0, ns, ts, ws, optODE );
 ceq = []; 
 cd = (0.5/f(3))*((f(1)+cb_in-f(2))*f(3) - (cao+cb_in-cbo)*Vo);
 
 c(1) = f(2) - cb_max;
 c(2) = cd - cd_max;
 else              % caso tenha análise de sensibilidade...
 
 [f,df] = fun( x0, ns, ts, ws, optODE );
 
 cd = (0.5/f(3))*((f(1)+cb_in-f(2))*f(3) - (cao+cb_in-cbo)*Vo);
 
 c(1) = f(2) - cb_max;
 c(2) = cd - cd_max;
 dcd = (1/2)*df(1,:) + (cao+cb_in-cbo)*Vo*df(3,:)/(f(3)^2)  - (1/2)*df(2,:);
 dcd = dcd';
 dc = df(2,:)';
 dc = [dc dcd];
 ceq = [];  
 dceq = [];  
   
end
end