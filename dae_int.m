function [t_s z] = dae_int(solver,function_handle,dt, z0, optODE )
    % This function integrates da dae model according to a specific solver
            [tspan zs] = ode45( function_handle,[dt], [z0], optODE );
            t_s = tspan;
            z = zs;
end
