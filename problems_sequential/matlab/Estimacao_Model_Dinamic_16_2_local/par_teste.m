clear all

contador=0;

for i=1:4153
    try
    load([num2str(i),'.mat'])
        if  J<420
            contador=contador +1;
            par_t{contador}=par;
        end
    catch ME
        J=1e30;
        par=[];
        save([num2str(i),'.mat'],'J','par')
        
    end
end

save('par_teste.mat','par_t')