function modelo = constraints(x,modelo,km_glyc,km_lact,vmax_glyc,vmax_lact,v_max_o2,km_o2,v_lactato)

% modelo.c = zeros(size(modelo.c));
% objectiveCoeff=1;
% 
% modelo.c(8)=objectiveCoeff;

if x(3)>0.001

    modelo = changeRxnBounds(modelo,'EX_glyc(e)', -vmax_glyc*(x(3)/(x(3)+km_glyc)),'l');

%    modelo = changeRxnBounds(modelo,'EX_ac(e)', 0,'l');

%    modelo = changeRxnBounds(modelo,'EX_o2(e)', -20*(x(9)/(x(9)+0.024)),'l');
	
%    modelo = changeRxnBounds(modelo,'EX_lac_D(e)', 0,'l');
else

%    modelo = changeRxnBounds(modelo,'EX_ac(e)',-2,'l');
    
    modelo = changeRxnBounds(modelo,'EX_glyc(e)',0,'l');
 
 %   modelo = changeRxnBounds(modelo,'EX_lac_D(e)', -2*(x(5)/(x(5)+0.01)),'l');
 
%    if x(4)>0.01
        
%    modelo = changeRxnBounds(modelo,'EX_ac(e)',-2*(x(4)/(x(4)+0.01)),'l');
	 
 %   else
 %        modelo = changeRxnBounds(modelo,'EX_ac(e)', 0,'l');
     
  %  end    
end

if x(9) > 0.0000001
    modelo = changeRxnBounds(modelo,'EX_o2(e)', -v_max_o2*(x(9)/(x(9)+km_o2)),'l');
%     modelo = changeRxnBounds(modelo,'EX_ac(e)',0,'l');
else
    modelo = changeRxnBounds(modelo,'EX_o2(e)', 0,'l');
%     modelo = changeRxnBounds(modelo,'EX_ac(e)',0,'l');
end

if x(6)>0.0001
    modelo = changeRxnBounds(modelo,'EX_lcts(e)',-vmax_lact*(x(6)/(x(6)+km_lact)),'l');
    modelo = changeRxnBounds(modelo,'EX_lac-D(e)',v_lactato,'l');
    modelo = changeRxnBounds(modelo,'EX_ac(e)',0,'l');
else
    modelo = changeRxnBounds(modelo,'EX_lcts(e)',-0.1,'l');
    modelo = changeRxnBounds(modelo,'plasm_lactose',0,'b');
    modelo = changeRxnBounds(modelo,'Ex_plasmideo_ast',0,'b');
    modelo = changeRxnBounds(modelo,'Replication_plasm',1e-6,'b'); 
    modelo = addExchangeRxn(modelo,{'plasmideo'},-20,20);
    modelo = changeRxnBounds(modelo,'EX_lac-D(e)',0,'b');
    modelo = changeRxnBounds(modelo,'EX_ac(e)',0,'b');
end

    