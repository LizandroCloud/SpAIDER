

 wname = 'db4';  % tipo de base wavelets.
 limiar = 1;
 w = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 10 10 10 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %sinal original
 w = w +rand(1,length(w))*0.8; % adiciona um pouco de ruído
 maxlev = wmaxlev(length(w),wname); % nível 
 [C,L] = wavedec(w,maxlev,wname);  % C sao os "pacotes de" detalhes e L  dimensão de cada C...
 % Obs: o primeiro C é aproximação, os outros são detalhes.
 for i=L(2)+1:length(C)
 if abs(C(i))<0.2*norm(C(2:end));  %limiar heuristico
    C(i)=0.0; % zera coeficiente
 end      
 end
 w_back = waverec(C,L,wname); % transformada wavelets  
 plot(w,'r'); hold on; plot(w_back, 'b')
                          
           
 