function [r_ix4b] =  g_ix4b(ix4,ic,C,is,nUt,difd,ns_a)
    figure(ix4)
    subplot(2,1,2) 
%  hist(difd(ic,1:length(C),is-1));   
 hist(difd);
 if nUt==1
 title('Distribution of details smaller than threshold'); 
 else
     title(['Distribution of details smaller than threshold  - control variable: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
 end
 ylabel('Values of details')
 xlabel('Frequence')
 r_ix4b=1;
end
