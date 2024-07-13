function dxdt = ode_with_mccormick(t, x)
    % ODE system with bilinear term
    a = 1; b = 2; c = 0.5; % Example coefficients
    u = x(1) * x(2); % Bilinear term
    v = x(1) + x(2) - 1; % Another term
    
    % McCormick envelopes for the bilinear term
    u_lower = min(x(1)^2, x(2)^2) * min(x(1), x(2));
    u_upper = max(x(1)^2, x(2)^2) * max(x(1), x(2));

    u_lower = u;
    u_upper = u;
    
    % Compute ODEs with McCormick envelopes
    dxdt(1,1) = -a * x(1) + u_lower; % Lower bound
    dxdt(2,1) = -b * x(2) + c * u_upper + v; % Upper bound
end

