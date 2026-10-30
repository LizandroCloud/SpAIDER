load('iJO1366_PM_asp_plasm_novo_1');
modelo=iJO1366_PM_asp_plasm_novo_1;
[v,obj,erro]=fba2(modelo,[],[],50,1);
obj
ID_rxn_lactose = findRxnIDs(modelo,'EX_lcts(e)');
modelo.lb(ID_rxn_lactose)
modelo.lb(ID_rxn_lactose)=0;
[v,obj,erro]=fba2(modelo,[],[],50,1);
obj
[v,obj,erro]=fba2_multi(modelo,0.8,obj,50,1e-6,2)
obj
