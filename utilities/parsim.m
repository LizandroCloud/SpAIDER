function [t,x,u] = parsim(model,tspan,time, u )
    [t,x,u]= sim(model,tspan,[],[time cell2mat(u) ]);
    
end