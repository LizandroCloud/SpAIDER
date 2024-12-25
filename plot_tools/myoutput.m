function stop = myoutput(x,optimvalues,state)
        stop = false;
        if state == 'iter'
%           history = [history; x];
%          figure(199)
%          bar(optimvalues.iteration,abs(optimvalues.fval));
        hold on;
        elseif state == 'init'
          history = 0;
        end
        

        
    end