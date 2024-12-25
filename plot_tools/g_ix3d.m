    function [r_ix3d] = g_ix3d(ix3,is,is0,ic,thrs,L,C)
         figure(ix3);
         subplot(2,1,2)
         ko=1;
         barcolor = 'w';
         for i=2:(length(L)-1)
            ts_x1=(ko+1:ko+L(i));
            if (i>1)
               theshold=ones(1,L(i))*thrs(ic,i-1,is-is0);
               bar(ts_x1',theshold,0.5,barcolor);
               bar(ts_x1',-(theshold),0.5,barcolor);
               ko = ko+L(i);
            end 
            hold on;
         end
%          axis([0 length(C)+1 -max(C)*1.05  max(C)*1.05 ])
         r_ix3d=1;
    end