function [ID,f]=fluxos(modelo,v)

% Oxigênio

ID_o2 = findMetIDs(modelo,'o2[e]'); 

% ID_rxn_o2 = findRxnIDs(modelo,'EX_o2(e)');

ID_rxn_o2 = findRxnIDs(modelo,'O2tex');

v_o2 = -v(ID_rxn_o2)*modelo.S(ID_o2,ID_rxn_o2);

%Substrato: Glicerol

ID_glyc = findMetIDs(modelo,'glyc[e]') ;

ID_rxn_glyc = findRxnIDs(modelo,'EX_glyc(e)');

v_glyc = -v(ID_rxn_glyc)*modelo.S(ID_glyc,ID_rxn_glyc);

%Indutor: Lactose
ID_lactose = findMetIDs(modelo,'lcts[e]') ;

ID_rxn_lactose = findRxnIDs(modelo,'EX_lcts(e)');

v_lactose = -v(ID_rxn_lactose)*modelo.S(ID_lactose,ID_rxn_lactose);


% Produtos extracelulares
    %Acetato
    
    ID_ac = findMetIDs(modelo,'ac[e]') ;

    ID_rxn_ac = findRxnIDs(modelo,'EX_ac(e)');
    
    v_ac = -v(ID_rxn_ac)*modelo.S(ID_ac,ID_rxn_ac);

    
    %Lactato
    
    ID_lac = findMetIDs(modelo,'lac-D[e]');

    ID_rxn_lac = findRxnIDs(modelo,'EX_lac-D(e)');
    
    v_lac = -v(ID_rxn_lac)*modelo.S(ID_lac,ID_rxn_lac);
    
    % Asparaginase
    
    ID_asp = findMetIDs(modelo,'asp');

    ID_rxn_asp = findRxnIDs(modelo,'asparaginase');
    
    v_asp = v(ID_rxn_asp)*modelo.S(ID_asp,ID_rxn_asp);
    
    % Plasmídeo
    
    ID_plasm = findMetIDs(modelo,'plasmideo');

    ID_rxn_plasm = findRxnIDs(modelo,'Replication_plasm');
    
    v_plasm = v(ID_rxn_plasm)*modelo.S(ID_plasm,ID_rxn_plasm);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ID=[ID_glyc ID_ac ID_lac ID_lactose ID_asp ID_plasm ID_o2];
    f=[v_glyc v_ac v_lac v_lactose v_asp v_plasm v_o2];
