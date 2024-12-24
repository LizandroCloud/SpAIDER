function [Pgrad] = adapt(w01,C,L,maxlev,ic,spp,ts,adap,normfrac)

    Pgrad=0;
    Pgrad = Blips(w01(ic,:),C,L,maxlev) % onlly selection of refined locations

    if isempty(Pgrad)
        Pgrad=0;
    end
    
Is=1;
while Is==1,
% for i=0:Pgrad(2)
if Pgrad(end) > length(w01(ic,:))
    Pgrad = Pgrad(1:end-1);
    Is = 1;
else
    Is = 0;
end
% end
% Is = 1;
end

if ( (Pgrad)==0 ) 
    return
end

if ( isempty(Pgrad)) 
    return
end
    dtA_i_max = 0; dtA_s_max = 0;
    [PgradN dtA] = Pgrad_ref(Pgrad,spp,ts,ic);
    dtA_max = max(dtA);
    
    if PgradN==0
        Pgradi = Pgrad -1;
        [PgradN dtA_i]= Pgrad_ref(Pgradi,spp,ts,ic);
        dtA_i_max = max(dtA_i);
    end
    
    
    if PgradN == 0
        Pgrads = Pgrad + 1;
        [PgradN dtA_s]= Pgrad_ref(Pgrads,spp,ts,ic);
        dtA_s_max = max(dtA_s);
    end
    
   
        dtA_g_max = max([dtA_max dtA_i_max dtA_s_max]);
        if dtA_g_max == dtA_i_max 
            kcl = -1;
        elseif dtA_g_max == dtA_max
            kcl = 0;
        elseif  dtA_g_max == dtA_s_max 
            kcl = 1;
        end
    Pgrad = Pgrad + kcl;
    spp=(dtA_g_max)^(-1);
    while PgradN==0
        spp=spp*(1.11);
        [PgradN dtA] = Pgrad_ref(Pgrad,spp,ts,ic);
    end
    
Pgrad = PgradN;
end