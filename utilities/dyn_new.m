function [param eq eq_sens n_p n_c n_a y J dJ cr_eq d_cr_eq cr_deq d_cr_deq cr_eq_type cr_deq_type n_s] = dyn 
 % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www.ceopt.com
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 %                              Description
 % Routine to calculate Jacobian of DAE model, objective function and
 % constraints
 %
 %*************************************************************************
 % changes:
 % 01-18-2012 - Calculation of Obj function and constraints jacobians
 %
 %*************************************************************************
 % For more information, please read the article:
 % http://www.sciencedirect.com/science/article/pii/B9780444595201501032
 % ************************************************************************
 dA='';  %pre-allocating
[ dx y par e_var a_var c_var J cr_eq cr_deq p] = input_model; % calling dae_model

[e_var a_var c_var] =  dae_var(p.nU,dx,y);
%    
param=char(par); % vector of (strings) parameters to be evaluated in model function
eq = char(dx);  % vector of (strings) equations to be evaluated in model function
eqo = eq;

    % symbolic
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25
    syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 tout
     syms y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16 y17 y18 y19 y20 y21 y22 y23 y24 y25 % no use **

% Sensibility equations
%**************************************************************************
if isempty(char(a_var))  % in case of only ode
    [R1] = jacmat_prior( eq , e_var, par, c_var ); % derivatives in relation to state var and control var
    eq_sens = char(R1); % vector of (strings) sensibility equations to be evaluated in model function
    n_p = length(par); % number of parameters
    n_c = length(c_var); % number of controls
    n_a = 0;
    y = char(y);
else % in case of DAE
    for j=1:length(e_var)
    for k=1:5
    for i=1:length(a_var) % re-arranging (implicit derivation)
        s1 = regexptranslate('escape',char(a_var(i))); %translate to calculate...
        a_eq  = regexprep((eqo(1,:)), s1,  (y(i)) ); % y() instead of y
        eqo = a_eq; % storage
    end
    end
    eq1(j,1:length(eqo)) = a_eq; % storage
    eqo=[];  eqo=eq(j,:);
    end
    [R2] = jacmat_prior( char(eq1) , e_var, par, c_var ); % derivatives in relation to state var and control var
    eq_sens = char(R2); % vector of (strings) sensibility equations to be evaluated in model function
    n_p = length(par); % number of parameters
    n_c = length(c_var); % number of controls
    n_a = length(a_var); % number of agebraic equations
    y = char(y);
end
% 
% Objective function gradient
%*************************************************************************
dJ=''; %preallocating
J1  = regexprep(J, 'end,', ''); % treating for gradient calculation...
for i=1:length(e_var)  % calculate for each state variable
    a=int2str(i); % converting to int
    [R3] = jacmat_states( char(J1) , e_var(i), par, c_var ); % derivatives in relation to state var
  %  R3 = regexprep(R3, 'x', 'f');  % change x for f to print later
	        for j=length(e_var):-1:1 % path contraint 
            b=int2str(j);
            sk = regexptranslate('escape',char(e_var(j))); %translate to calculate...
            R3 = regexprep(R3, sk, ['f(end,',b,')']); % changing
            end
    dA=['(',char(R3),')','*df(end,',a,',:)']; % substituting again to print later
    if str2num(a)=='1' 
        dJ = [dA]; % first iteration
    else
        dJ = [dJ,'+',dA]; % other iterations
    end
end
for j=1:length(c_var)
[R4] = jacmat_control( char(J1) , e_var(i), par, c_var(j) ); % derivatives in relation to control var
dJ = [dJ,'+',char(R4),'*ones(1,1,length(ws))']; % adding the control derivative part
end
%dJ = [dJ,'+',char(R4),'*ones(1,1,length(ws))']; % adding the control derivative part
J = regexprep(J, 'x', 'f');  % change x for f to print later
J=char(J); % converting to char

% Equality constraints...
%*************************************************************************
if length(char(cr_eq))~=2

cr_eq_type=ones(1,length(cr_eq))*0; % pre-allocating
d_cr_eq=blanks(100); % pre-allocating
for j=1:length(cr_eq)
    cr_eq_1  = regexprep(cr_eq(j), 'end,', ''); % treating for gradient calculation...
    b=strcmp(cr_eq_1,char(cr_eq(j))); % comparing strings
    if b==0 % if exists the word 'end';
        cr_t = 'final';
    else
        cr_t = 'path';
    end
    cr_eq_2  = regexprep(cr_eq_1, ':,', ''); % treating for gradient calculation...
%     [R5] = jacmat_prior( char(cr_eq_2) , e_var, par, c_var ); 
    for i=1:length(e_var) % calculate for each state variable
        a=int2str(i);  % converting to int
        [R5] = jacmat_states( char(cr_eq_2) , e_var(i), par, c_var );  % derivatives in relation to states var
        cr_eq_type(j)=strcmp(cr_t,'final'); % comparing to 'final'
            if cr_eq_type(j)==1 % if cr_type == 'final'
                dA=['(',char(R5),')','*df(end,',a,',:)']; % putting 'end'
                    if a=='1' % for the first iteration
                        d_cr_eq = dA; % writing
                    else % other iterations
                        d_cr_eq = [d_cr_eq,'+',dA]; % writing
                    end
            else  % if cr_type ~= 'final'
                dA=['(',char(R5),')','*df(time,',a,',:)'];
                    if a=='1'% for the first iteration
                        d_cr_eq = [dA]; % writing
                    else % other iterations
                        d_cr_eq = [d_cr_eq,'+',dA]; % writing
                    end                
            end
    end
    dB(j,:) = {d_cr_eq}; % writing
    if cr_eq_type(j)==1 % tf contraint 
        for i=1:length(e_var)
            a=int2str(i);
            sk = regexptranslate('escape',char(e_var(i))); %translate to calculate...
            dB(j) = regexprep(dB(j), sk, ['f(end,',a,')']); % changing
        end
    else
        for i=1:length(e_var) % path contraint 
            a=int2str(i);
            sk = regexptranslate('escape',char(e_var(i))); %translate to calculate...
            dB(j) = regexprep(dB(j), sk, ['f(time,',a,')']); % changing
        end
    end
    cr_eq(j) = regexprep(cr_eq(j), 'x(', 'f('); % changing 
    cr_eq(j) = regexprep(cr_eq(j), ':', 'time'); % changing 
end
cr_eq = char(cr_eq); % converting to string
d_cr_eq = char(dB);   % converting to string
else
cr_eq = char({'[]'}); % converting to string
d_cr_eq = char({'[]'});   % converting to string     
cr_eq_type=99;  % empty

end


% Inequality constraints...
%*************************************************************************
if length(char(cr_deq))~=2
cr_deq_type=ones(1,length(cr_deq))*0; % pre-allocating
d_cr_deq=blanks(100); % pre-allocating
for j=1:length(cr_deq)
    cr_deq_1  = regexprep(cr_deq(j), 'end,', ''); % treating for gradient calculation...
    b=strcmp(cr_deq_1,char(cr_deq(j))); % comparing strings
    if b==0 % if exists the word 'end';
        cr_t = 'final';
    else
        cr_t = 'path';
    end
    cr_deq_2  = regexprep(cr_deq_1, ':,', ''); % treating for gradient calculation...
%     [R5] = jacmat_prior( char(cr_deq_2) , e_var, par, c_var ); 
    for i=1:length(e_var) % calculate for each state variable
        a=int2str(i);  % converting to int
        [R6] = jacmat_states( char(cr_deq_2) , e_var(i), par, c_var );  % derivatives in relation to states var
        cr_deq_type(j)=strcmp(cr_t,'final'); % comparing to 'final'
            if cr_deq_type(j)==1 % if cr_type == 'final'
                dA=['(',char(R6),')','*df(end,',a,',:)']; % putting 'end'
                    if a=='1' % for the first iteration
                        d_cr_deq = dA; % writing
                    else % other iterations
                        d_cr_deq = [d_cr_deq,'+',dA]; % writing
                    end
            else  % if cr_type ~= 'final'
                dA=['(',char(R6),')','*df(time,',a,',:)'];
                    if a=='1'% for the first iteration
                        d_cr_deq = [dA]; % writing
                    else % other iterations
                        d_cr_deq = [d_cr_deq,'+',dA]; % writing
                    end                
            end
    end
    dC(j,:) = {d_cr_deq}; % writing
    if cr_deq_type(j)==1 % tf contraint 
        for i=1:length(e_var)
            a=int2str(i);
            sk = regexptranslate('escape',char(e_var(i))); %translate to calculate...
            dC(j) = regexprep(dC(j), sk, ['f(end,',a,')']); % changing
        end
        dC(j) = regexprep(dC(j), 'y(', 'y(end,'); % changing
    else
        for i=1:length(e_var) % path contraint 
            a=int2str(i);
            sk = regexptranslate('escape',char(e_var(i))); %translate to calculate...
            dC(j) = regexprep(dC(j), sk, ['f(time,',a,')']); % changing
            dC(j) = regexprep(dC(j), 'x(', 'f('); % changing
        end
        dC(j) = regexprep(dC(j), 'y(', 'y(time,'); % changing
    end
    cr_deq(j) = regexprep(cr_deq(j), 'x(', 'f('); % changing 
    cr_deq(j) = regexprep(cr_deq(j), ':', 'time'); % changing 
   
end
    cr_deq = char(cr_deq); % converting to string
    d_cr_deq = char(dC);   % converting to string
    else
cr_deq = char({'[]'}); % converting to string
d_cr_deq = char({'[]'});   % converting to string    
cr_deq_type=99;  % empty

end

n_s = length(e_var);
end
