function h = haar2(n,j,k,t)

if t==1

t = (2^j)*t - k;
% Escrevendo funcao haar...
%++++++++++++++++++++++++++

if t>=0
    if t<0.5
    h=1;  % ok...
    else  % se for maior q 0.5..
        if t<=1  % mas menor ou igual a 1...
            h=-1;  %ok...
        else  % se passar de 1
            h=0;
        end
    end
else
    h=0.0;  
end



else
    
t = (2^j)*t - k;
% Escrevendo funcao haar...
%++++++++++++++++++++++++++

if t>=0
    if t<0.5
    h=1;  % ok...
    else  % se for maior q 0.5..
        if t<1  % mas menor ou igual a 1...
            h=-1;  %ok...
        else  % se passar de 1
            h=0;
        end
    end
else
    h=0.0;  
end

end
end
    