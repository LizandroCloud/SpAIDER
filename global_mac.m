function [a] = global_mac()
    % Define time span
    tspan = [0 10];
    
    % Initial condition
    x0 = [1; 0.5]; % Example initial condition
    
    % Solve ODEs with McCormick envelopes
    [t, x] = ode45(@ode_with_mccormick, tspan, x0);
    
    % Plot results
    plot(t, x, '--');
    xlabel('Time');
    ylabel('State');
    legend('x1', 'x2');
    title('ODEs with McCormick Envelopes');
    grid on;
    a=0;
end
