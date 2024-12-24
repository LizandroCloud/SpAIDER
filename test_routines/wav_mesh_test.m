function [cfd,detais0,maxlev,C,L] = wav_mesh_test(wopt,is,ic,wname)
 % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: 
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 %                              Description
 % Function to generate wavelets mesh...
 %
 % 
 % ************************************************************************
 
 % changes:
 % 01-21-2013 - First modification
 %*************************************************************************
maxlev = wmaxlev(length(wopt),wname); % maximum level of wavelet resolution 
%    wopt1 = wopt*2^(-maxlev/2);
[C,L] = wavedec(wopt,maxlev,wname);  % C are wavelets coefficients
detais0(ic,1:length(C),is)=C;  % details before thresholding

for i=1:length(C)
    ci(1,i)=C(i);
end

for i=1:length(L)
    li(1,i)=L(i);
end

cfd = zeros(maxlev,length(wopt));
for k = 1:maxlev % extracting wavelets details for plotting
     d = detcoef(ci,li,k); 
     d = d(ones(1,2^k),:);   
    cfd(k,:) = wkeep(d(:)',length(wopt)); 
end 

cfd = cfd(:); 
I = find(abs(cfd)<sqrt(eps)); 
cfd(I)=zeros(size(I)); 
cfd = reshape(cfd,maxlev,length(wopt));

end
