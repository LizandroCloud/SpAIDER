
% select the cvodesnonstiff solver
d = ode(ODEFcn=@(t,y,p) [p(1)*y(1) - p(2)*y(1)*y(2); ...
                           p(3)*y(1)*y(2) - p(4)*y(2)]);

% Set initial conditions and parameters
d.InitialTime = 0;
d.InitialValue = [1;1];
d.Parameters = [.4 .4 .09 .1];
 
% Set the relative tolerance
d.RelativeTolerance = 1e-6;

% Select the ode45 solver
d.Solver = "auto";
 
% Solve over a short time span
sol = solve(d,0,100);
plot(sol.Time,sol.Solution,LineWidth=2)