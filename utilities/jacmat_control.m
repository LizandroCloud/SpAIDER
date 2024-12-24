function [j4] = jacmat_control(dfdu , e_var , par, c_var )

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

for i=1:length(dfdu(:,1))
  
  for j=1:length(c_var)
      if j==1
            s1 = regexptranslate('escape',char(c_var(j)));
            s2 = regexprep(dfdu(i,:), s1, [ ' u',int2str(j),' ']);
            dfdu1(i,1:length(s2)) = s2;
      else
            s1 = regexptranslate('escape',char(c_var(j)));
            s2 = regexprep(dfdu1(i,:), s1, [ ' u',int2str(j),' ']);
            dfdu1(i,1:length(s2)) = s2;
      end
  end
  
 end 


 
 for i=1:length(dfdu1(:,1))
  
  for j=1:length(e_var)
      if j==1
            s1 = regexptranslate('escape',char(e_var(j)));
            s2 = regexprep(dfdu1(i,:), s1, [ ' x',int2str(j),' ']);
            dfdu2(i,1:length(s2)) = s2;
      else
            s1 = regexptranslate('escape',char(e_var(j)));
            s2 = regexprep(dfdu2(i,:), s1, [ ' x',int2str(j),' ']);
            dfdu2(i,1:length(s2)) = s2;
      end
  end
   dfdu10(i,1:length(s2)) = sym(str2sym(dfdu2(i,1:length(s2))));
   s2=0;
 end 

 for j=1:length(e_var)
      
            s1 = regexptranslate('escape',char(e_var(j)));
            s2 = regexprep(e_var(j), s1, [ ' x',int2str(j),' ']);
            e_var1(j) = sym(char(s2));
 end
 
 for j=1:length(c_var)
      
            s1 = regexptranslate('escape',char(c_var(j)));
            s2 = regexprep(c_var(j), s1, [ ' u',int2str(j),' ']);
            c_var1(j) = sym(char(s2));
 end
 
 
for i=1:length(dfdu1(:,1))
    df(i,:)=jacobian(str2sym(dfdu2(i,:)),[c_var1]);  % not jacobian2
end


% for j=1:length(dfdu1(:,1))
%     for i=1:length(v)
%         dfduu(j,:)=subs(df( j,: ),v(i),x(i),0);
%         df(j,:) = dfduu(j,:);
%     end
% % jac(j) = sum(df(j,:))';
% end
%  jac = (df);
% jac=roundn(jac,-2);
% % j2 = char(ones(100,100));

for i=1:length(dfdu(:,1))
      s=char(df(i,:));
      df1 = regexprep(s, 'matrix', '');
  for j=1:length(c_var)

      if j==1
            s1 = regexptranslate('escape',char(c_var(j)));
            s2 = regexprep(df1, [ 'u',int2str(j)], s1);
            df2(i,1:length(s2)) = s2;
      else
            s1 = regexptranslate('escape',char(c_var(j)));
            s2 = regexprep(df2(i,1:length(s2)), [ 'u',int2str(j)], s1);
            df2(i,1:length(s2)) = s2;
      end
  end
  df3(i,:) = {df2(i,:)};
end


for i=1:length(dfdu(:,1))
      df4=char(df3(i,:));
%       df1 = regexprep(s, 'matrix', '');
  for j=1:length(e_var)

      if j==1
            s1 = regexptranslate('escape',char(e_var(j)));
            s2 = regexprep(df4, [ 'x',int2str(j)], s1);
            df5(i,1:length(s2)) = s2;
      else
            s1 = regexptranslate('escape',char(e_var(j)));
            s2 = regexprep(df5(i,1:length(s2)), [ 'x',int2str(j)], s1);
            df5(i,1:length(s2)) = s2;
      end
  end
  df6(i,:) = {df5(i,:)};
end
  j2=df6;
  j3  = regexprep(df6, '([', '');
  j4  = regexprep(j3, '])', '');

end