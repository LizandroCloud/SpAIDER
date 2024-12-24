function h = haar(t,J)

h=[];
h(1)=1;  % funcao fundamental...

% Escrevendo funcao haar...
%++++++++++++++++++++++++++
% if t==1
%     h(1)=1;
% end
% 
% if t==0
%     h(1)=1;
% end

for j = 0:J  % funcao para os outros niveis...
    for k = 0:(2^j - 1);
        n = k + 2^j + 1;
        hi=haar2(n,j,k,t);
        h(n)=hi;
        hi=0.0;
    end
end
