function r_ix1 = g_ix1(ix1,ic,ns,ns_a,cfd,nUt,maxlev)
figure(ix1)
subplot(2,1,1), colormap(1-gray(128)); 
imagesc(wcodemat(flipud(cfd(1:maxlev,1:ns(ic)))));
% grid on;
% img = imagesc((wcodemat(cfd,128,'r'))); 
% set(get(img,'parent'),'YtickLabel',[]); 

if nUt==1
    title(['Original wavelets coefficients for: ',int2str(ns_a(ic,:)),' stages']); 
else
    title(['Original wavelets coefficients - control variable: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
end
ylabel('Resolution levels')
xlabel('Mesh length - Control Stages')
r_ix1=1;
end