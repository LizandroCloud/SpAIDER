function pfd = DecryptSystem(y,vars,eqps,pfd)
    %Try later build a class for pfd
    for k=1:size(y,2)
        for j=1:vars.NumOfVars
            if k == vars.Variables{j}.ID
                evalc(strcat('pfd.',vars.Variables{j}.Name,' = vars.Variables{j}'));
                evalc(strcat('pfd.',vars.Variables{j}.Name,'.Value = y(:,k)'));
            end
        end
        for j=1:eqps.NumOfEqps
            for m=1:eqps.Equipment{j}.NumOfVars
                if k == eqps.Equipment{j}.Variables{m}.ID
                    evalc(strcat('pfd.',eqps.Equipment{j}.Variables{m}.Name,' = eqps.Equipment{j}.Variables{m}'));
                    evalc(strcat('pfd.',eqps.Equipment{j}.Variables{m}.Name,'.Value = y(:,k)'));
                    %Try later to use the struct command
                end
            end
        end
    end 
end

