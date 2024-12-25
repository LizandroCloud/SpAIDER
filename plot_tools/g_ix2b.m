function [r_ix2b] = g_ix2b(ix2,ic,ns_a,talf,C,nUt,L)
C_x=1;
figure(ix2);
C_x(1:length(C)) = C(1:length(C));
ko=0;
subplot(2,1,2);
tss = [1:length(C_x)];
axis([0 length(C) 1.2*min(C)  max(C)*1.2 ]);
stairs(tss*talf,zeros(1,length(C_x)),'w');
% legend('Ap')
hold on;
for i=1:(length(L)-1)  
     figure(ix2);
     subplot(2,1,2); 
      axis([0 length(C)+1 1.2*min(C)  max(C)*1.2 ]) 
    if i==1 
        barcolor = 'k'; 
        elseif i==2 barcolor = 'b'; 
        legend('Ap','Dt.1')
        elseif i==3 barcolor = 'g';
        legend('Ap','Dt.1','Dt.2')
        elseif i==4 barcolor = 'r'; 
        legend('Ap','Dt.1','Dt.2','Dt.3')
        elseif i==5 barcolor = 'y';
        legend('Ap','Dt.1','Dt.2','Dt.3','Dt.4')
        elseif i==6 barcolor = 'm';   
        legend('Ap','Dt.1','Dt.2','Dt.3','Dt.4','Dt.5')
        elseif i==7 barcolor = 'c';   
    end
    
    ts_x1=(ko+1:ko+L(i));
    C_x1 = ( C_x(ko+1:ko+L(i)) );
    if length(ts_x1) == 1
        ts_x1 = [ts_x1 ( ko+ L(i)+1)];
        C_x1 = [(C_x1) (C_x( ko+ L(i)+1))];
    end
    
    bar(ts_x1,C_x1,0.8,'stacked',barcolor); 
    axis([0 length(C)+1 1.2*min(C)  max(C)*1.2 ])
    ko = ko+L(i);
    hold on;
end
grid on;
i=length(L)-1;
if i==1 
        elseif i==2 
        legend('Legend','Ap','Dt.1')
        elseif i==3 
        legend('Legend','Ap','Dt.1','Dt.2')
        elseif i==4 
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3')
        elseif i==5
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3','Dt.4')
        elseif i==6 
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3','Dt.4','Dt.5')
        elseif i==7 
        legend('Legend','Ap','Dt.1','Dt.2','Dt.3','Dt.4','Dt.5','Dt.6')
    end
if nUt==1
 title('Wavelets aproximation and details coefficients'); 
 else
     title(['Wavelets aproximation and details coefficients  - control variable: ',int2str(ic),' for: ',int2str(ns_a(ic,:)),' stages'])
 end
ylabel('Wavelets coefficients')
xlabel('Points')
r_ix2b=1;
end