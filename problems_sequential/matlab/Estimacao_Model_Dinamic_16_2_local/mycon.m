 function [c,ceq] = mycon(v,mi_ast)
         c = [];
         %Considerando alfa=0.6
         ceq = v(13)-0.6*mi_ast;
         
 