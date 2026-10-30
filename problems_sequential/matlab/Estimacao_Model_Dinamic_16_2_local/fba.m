function [v,obj]=fba(modelo,num)

% modelo = changeRxnBounds(modelo,'Ec_biomass_iJO1366_WT_53p95M',0,'b');
% modelo = changeRxnBounds(modelo,'Ec_biomass_iAF1260_core_59p81M',0,'b');
modelo = changeRxnBounds(modelo,'EX_glc(e)',0,'b');
modelo = changeRxnBounds(modelo,'EX_glyc(e)',(-8),'l'); 
modelo = changeRxnBounds(modelo,'NADTRHD',4,'u'); % NAD transhydrogenase
modelo = changeRxnBounds(modelo,'G6PDH2r',3.5,'u'); % glucose 6 phosphate dehydrogenase - via das pentoses

modelo = changeObjective(modelo,{'Ec_biomass_iJO1366_core_53p95M'});
% modelo = changeObjective(modelo,{'Ec_biomass_iAF1260_core_59p81M'});

if num==1
% FBA clássico
changeCobraSolver('glpk');

FBAsolution = optimizeCbModel(modelo,'max');

v=FBAsolution.x;
obj = FBAsolution.f;

switch FBAsolution.stat
    case  1 
        disp('Optimal solution')
    case 2 
        disp('Unbounded solution')
    case 0
        disp('Infeasible')
    case -1
        disp('No solution reported (timelimit, numerical problem etc)')
end

%FBA não-clássico
elseif num==2
    
modelo = changeRxnBounds(modelo,'Ec_biomass_iJO1366_core_53p95M',0.2,'b');
% modelo = changeRxnBounds(modelo,'Ec_biomass_iAF1260_core_59p81M',0.2,'b');

changeCobraSolver('gurobi5');

c=zeros(length(modelo.c),1);

QPproblem= struct('A', modelo.S, 'lb', modelo.lb, 'ub', modelo.ub, 'b', zeros(size(modelo.S,1),1),'c',c,'osense', 0, 'csense', repmat('E',size(modelo.S,1),1),'F',(2)*diag(ones(length(modelo.c),1)));

solution = solveCobraQP(QPproblem,'default');

v=solution.full;
obj=solution.obj;

switch solution.stat
    case  1 
        disp('Optimal solution')
    case 2 
        disp('Unbounded solution')
    case 0
        disp('Infeasible')
    case -1
        disp('No solution reported (timelimit, numerical problem etc)')
end

elseif num==3
    op=optimset('Display','off');
    [v,FVAL,EXITFLAG] = linprog(-modelo.c',[],[],modelo.S,modelo.b,modelo.lb,modelo.ub,[],op);
    obj=-FVAL;
    
    switch EXITFLAG
          case 1
            disp('Linprog converged to a solution X')
          case 0
             disp('Maximum number of iterations reached.')
          case -2
             disp('No feasible point found.')
          case -3
             disp('Problem is unbounded.')
          case -4
             disp('NaN value encountered during execution of algorithm.')
          case -5
             disp('Both primal and dual problems are infeasible.')
          case -7
             disp('Magnitude of search direction became too small; no further progress can be made. The problem is ill-posed or badly conditioned.')
    end
%FBA não-clássico
elseif num==4
    options = optimset('Algorithm','active-set','Maxiter',1000,'display','off'); 
    H=diag(2*ones(length(modelo.c),1));
    Aeq=zeros(1,length(modelo.c));
    Aeq(8)=1;
    [i,j,s]=find(Aeq);
    [m,n]=size(Aeq);
    Aeq=sparse(i,j,s,m,n);
%     beq=0.6*mi_ast;
    beq=0.2;
    c=zeros(length(modelo.c),1);
    matrix=[modelo.S; Aeq];
    [v,fval,exitflag] = quadprog(H,c,[],[],full(matrix),[modelo.b; beq],modelo.lb,modelo.ub,[],options);
    obj=fval;
    switch exitflag
          case 1
            disp('Linprog converged to a solution X')
          case 0
             disp('Maximum number of iterations reached.')
          case -2
             disp('No feasible point found.')
          case -3
             disp('Problem is unbounded.')
          case -4
             disp('NaN value encountered during execution of algorithm.')
          case -5
             disp('Both primal and dual problems are infeasible.')
          case -7
             disp('Magnitude of search direction became too small; no further progress can be made. The problem is ill-posed or badly conditioned.')
    end
    
elseif num==5
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

elseif num==6
%FBA não-clássico
modelo = changeRxnBounds(modelo,'Ec_biomass_iJO1366_core_53p95M',0.2,'b');
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
   
    model.A = modelo.S;
    model.obj = zeros(length(modelo.c),1);
    model.sense = '=';
    model.rhs = modelo.b;
    
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
    model.constrnames = modelo.metNames;
    result=gurobi(model);
    v=result.x;
    obj=result.objval;
    
end