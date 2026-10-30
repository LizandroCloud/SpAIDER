function [jac,j1] = jacmat(dfdu , v , x )

 % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www. ????
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 %                              Description
 % Routine for calculation of sensitivity equations of dynamic model
 %  
 %
 %
 % ************************************************************************

% Variables and parameters: 
%
% [in]:
%
% dfdu --> diefferential equantions
% v --> independent variables
% x --> state
% [out]:
%
% j1  -->  sensitivity equations
% jac -->  sensitivity equations in corrent state x

for i=1:length(dfdu)
    df(i,:)=jacobian(dfdu(i),v);
end
j1 = df;
for j=1:length(dfdu)
    for i=1:length(v)
        dfduu(j,:)=subs(df( j,: ),v(i),x(i),0);
        df(j,:) = dfduu(j,:);
    end
% jac(j) = sum(df(j,:))';
end
jac = df;
jac=roundn(jac,-2);

end