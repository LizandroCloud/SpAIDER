function [e_var a_var c_var] =  dae_var(nU,x0,y)

 % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www.ceopt.com
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 %                              Description
 % Routine to write variables for jacobian calculation
 % [in]
 % nU  -> number of control variables
 % x0 -> initial state variables
 % y   -> algebraic equations
 % [out]
 % e_var -> state variables
 % a_var -> algebraic variables
 % c_var -> control variables
 %
 %*************************************************************************
 % changes:
 % 08-27-2013 - creation
 %
 %*************************************************************************
 % For more information, please read the article:
 % http://www.sciencedirect.com/science/article/pii/B9780444595201501032
 % ************************************************************************
     
        i=1:length(x0);
        letterx = repmat('x(',length(x0),1);
        e_vari = strcat(num2str(letterx),num2str(i'),')');
        e_varc = cellss (i,e_vari);
        e_var = regexprep(e_varc,'x( ','x(');
             
        j=1:nU;
        letteru = repmat('u(',nU,1);
        c_vari = strcat(num2str(letteru),num2str(j'),')');
        c_varc = cellss (j,c_vari);
        c_var = regexprep(c_varc,'u( ','u(');
        
        switch isempty(cell2mat(y))
            case{1}
                a_var={''};
            otherwise
                 k=1:length(y);
                 lettery = repmat('y(',length(y),1);
                 a_vari = strcat(num2str(lettery),num2str(k'),')');
                 a_varc = cellss (k,a_vari);
                 a_var = regexprep(a_varc,'y( ','y(');
        end
        
 
end
   