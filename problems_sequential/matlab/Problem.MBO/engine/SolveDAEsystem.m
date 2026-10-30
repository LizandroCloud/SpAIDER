%function [t,y] = SolveDAEsystem(resfile,init,vars,tspan,y0,ydot0,rtol,atol,solver)
function [t,y] = SolveDAEsystem(resfile,tspan,y0_guess,fixed_y0,ydot0_guess,fixed_ydot0,rtol,atol,solver)
    %CONSISTENT INITIAL CONDTIONS SOLVE
    tic
    disp(strjoin({'Starting constrained optimization to find consistent initial conditions for ',resfile,'...'},''));
    [y0,ydot0,rsnrm] = findcic(resfile,y0_guess,fixed_y0,ydot0_guess,fixed_ydot0,atol);
    elapsedTime = toc;
    if rsnrm < atol
        disp(strjoin({'Consistent initial conditions were obtained successfully in',num2str(elapsedTime),'seconds'}));
    else
        error(strjoin({'non-linear solver could not obtain consistent initial conditions. resnorm =',num2str(rsnrm)}));
    end
    
    %DAE SYSTEM SOLVE
    tic
    disp(strjoin({'Starting dynamic simulation of ',resfile,'...'},''));
    switch solver
        case 'ode15i'
            options = odeset(...
                'RelTol',rtol, ...
                'AbsTol',atol, ...
                'NormControl','on', ...
                'stats','on');
            [t,y] = ode15i(resfile,tspan,y0,ydot0,options);
        case 'dasslc'
            [t,y] = dasslc('resfile',tspan,y0,ydot0,rtol,atol);
    end
    elapsedTime = toc;
    disp(strjoin({'Simulation has finished successfully in',num2str(elapsedTime),'seconds'}));
end

