function [v,obj,erro]=fba1(modelo,k,mi,num)

if num==1
%     modelo = changeObjective(modelo,{'asparaginase'});
    modelo = changeObjective(modelo,{'Replication_plasm'});
%       modelo = changeObjective(modelo,{'Ec_biomass_iJO1366_core_53p95M'});
%     modelo = changeRxnBounds(modelo,'NADTRHD',0,'b'); % NAD transhydrogenase

%     modelo = changeRxnBounds(modelo,'EX_o2(e)',-5,'l'); % NAD transhydrogenase

%Resolução utilizando o solver gurobi5
%     minimize           c'*x + x'*Q*x + alpha
%        x
%     subject to         A*x = b,
%                        l <= x <= u,
%                        some xj integral,
%                        some xk must lie within second order cones,
%                        x'*Qc*x + q'*x <= beta,
%                        some xi in SOS constraints.
    
    model.A = modelo.S;
    model.obj = modelo.c;
    model.sense = '=';
    model.rhs = modelo.b;
    
    model.lb = modelo.lb;
    model.ub = modelo.ub;
    model.vtype = 'C'; %continuous variables 
    model.modelsense = 'max';
%     model.modelname = modelo;
%     model.objcon
%     model.vbasis
%     model.cbasis
%     model.Q
%     model.cones
%     model.quadcon
%     model.sos
%     model.pwlobj
%     model.start
    model.varnames = modelo.rxnNames;
    model.constrnames = modelo.metNames;
    result=gurobi(model);
    v=result.x;
    obj=result.objval;
    erro=max(modelo.S*v);
elseif num==2
%FBA não-clássico
modelo = changeRxnBounds(modelo,'Ec_biomass_iJO1366_core_53p95M',mi,'b');
% modelo = changeRxnBounds(modelo,'Ec_biomass_iAF1260_core_59p81M',0.2,'b');

        %Resolução utilizando o solver gurobi5
%     minimize           c'*x + x'*Q*x + alpha
%        x
%     subject to         A*x = b,
%                        l <= x <= u,
%                        some xj integral,
%                        some xk must lie within second order cones,
%                        x'*Qc*x + q'*x <= beta,
%                        some xi in SOS constraints.
    %Nova restrição
    %v_plasm=k*mi;
    Aeq=zeros(1,length(modelo.c));
    Aeq(2584)=1;
    model.A = [modelo.S;Aeq];
    model.obj = zeros(length(modelo.c),1);
    model.sense = '=';
    beq=k*mi;
    model.rhs = [modelo.b;beq];
    
    model.lb = modelo.lb;
    model.ub = modelo.ub;
    model.vtype = 'C'; %continuous variables 
    model.modelsense = 'min';
%     model.modelname = modelo;
%     model.objcon
%     model.vbasis
%     model.cbasis
    Q=diag(ones(length(modelo.c),1));
    [i,j,s]=find(Q);
    [m,n]=size(Q);
    model.Q = sparse(i,j,s,m,n); %Matriz tem que ser positiva semi-definida
%     model.cones
%     model.quadcon
%     model.sos
%     model.pwlobj
%     model.start
    model.varnames = modelo.rxnNames;
%     model.constrnames = modelo.metNames;
    result=gurobi(model);
    v=result.x;
    obj=result.objval;
    erro=unique(model.A*v);
    erro=erro(end-1);
end
