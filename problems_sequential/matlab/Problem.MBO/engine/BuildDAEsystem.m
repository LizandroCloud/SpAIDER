function [y0,y0_fixed,ydot0,ydot0_fixed] = BuildDAEsystem(daesys,vars,eqps,cnct,spec,eqs,init)
%DEFINING SYSTEM VARIABLES
    for k=1:eqps.NumOfEqps
        for j=1:eqps.Equipment{k}.NumOfVars
            AddVariable(vars,eqps.Equipment{k}.Variables{j}.Name,eqps.Equipment{k}.Variables{j});
        end
    end

%DEFINING SYSTEM VARIABLES IDENTIFIER
    idx = 0; y0 = [];
    for k=1:cnct.NumOfVarCnct
		id_VarFrom = GetVarIndex(vars,cnct.VarNameFrom{k});
		id_VarTo = GetVarIndex(vars,cnct.VarNameTo{k});
		if isequal(class(vars.Variables{id_VarFrom}),class(vars.Variables{id_VarTo})) == false;
			error('connection between different variable types is not possible');
		else
			if isempty(vars.Variables{id_VarFrom}.ID)
				idx=idx+1;
				vars.Variables{id_VarFrom}.ID = idx;
				y0 = [y0;vars.Variables{id_VarFrom}.Value];
			end
			vars.Variables{id_VarTo}.ID = vars.Variables{id_VarFrom}.ID;
		end
    end
    
    for k=1:cnct.NumOfEqpCnct
        EqpNameFromSplt = strsplit(cnct.EqpNameFrom{k},'.');
        EqpNameFrom = EqpNameFromSplt{1};
        id_EqpFrom = GetEqpIndex(eqps,EqpNameFrom);
        EqpStructFrom = strjoin({'eqps.Equipment{',num2str(id_EqpFrom),'}'},'');
        for j=2:length(EqpNameFromSplt)
            id_Eqp = GetEqpIndex(eval(EqpStructFrom),cnct.EqpNameFrom{k});
            EqpStructFrom = strjoin({EqpStructFrom,'.Equipment{',num2str(id_Eqp),'}'},'');
        end
        
        EqpNameToSplt = strsplit(cnct.EqpNameTo{k},'.');
        EqpNameTo = EqpNameToSplt{1};
        id_EqpTo = GetEqpIndex(eqps,EqpNameTo);
        EqpStructTo = strjoin({'eqps.Equipment{',num2str(id_EqpTo),'}'},'');
        for j=2:length(EqpNameToSplt)
            id_Eqp = GetEqpIndex(eval(EqpStructTo),cnct.EqpNameTo{k});
            EqpStructTo = strjoin({EqpStructTo,'.Equipment{',num2str(id_Eqp),'}'},'');
        end
        
        supcls = superclasses(eval(EqpStructFrom)); torf = true;
        for j=1:length(supcls)
            if (isequal(supcls{j},class(eval(EqpStructTo))) || isequal(class(eval(EqpStructFrom)),class(eval(EqpStructTo)))) && torf
                for m=1:eval(strjoin({EqpStructTo,'.NumOfVars'},''))
                    id_VarFrom = GetVarIndex(vars,eval(strjoin({EqpStructFrom,'.Variables{',num2str(m),'}.Name'},'')));

                    if isempty(vars.Variables{id_VarFrom}.ID)
                        idx=idx+1;
                        vars.Variables{id_VarFrom}.ID = idx;
                        y0 = [y0;vars.Variables{id_VarFrom}.Value];
                    end

                    id_VarTo = GetVarIndex(vars,eval(strjoin({EqpStructTo,'.Variables{',num2str(m),'}.Name',},'')));
                    if strcmp(vars.Variables{id_VarTo}.Name(length(cnct.EqpNameTo{k})+2:length(vars.Variables{id_VarTo}.Name)), ...
                            vars.Variables{id_VarFrom}.Name(length(cnct.EqpNameFrom{k})+2:length(vars.Variables{id_VarFrom}.Name)))
                        vars.Variables{id_VarTo}.ID = vars.Variables{id_VarFrom}.ID;
                    end
                end
                torf = false;
            end
        end
    end
	
	for k=1:vars.NumOfVars
		if isempty(vars.Variables{k}.ID)
            idx = idx+1;
            vars.Variables{k}.ID = idx;
            y0=[y0;vars.Variables{k}.Value];		
		end
	end

%DEFINING SYSTEM SIZE
    SysSize = idx; idx = 0;
    
%DEFINING SYSTEM EQUATIONS
    daeeqs = cell(SysSize,1);
    for k=1:spec.NumOfSpecs
        AddEquation(eqs,spec.Specifications{k})
    end
    
    for k=1:eqps.NumOfEqps
        for j=1:eqps.Equipment{k}.NumOfEqs
            AddEquation(eqs,eqps.Equipment{k}.Equations{j})
        end
    end
    
    for k=1:eqs.NumOfEqs
        idx = idx+1;
        fh = ReplaceVariable(eqs.Equations{k},vars);
        daeeqs{idx} = strjoin({'res(',num2str(idx),') = ',fh,';'},'');
    end    

    daefile = sprintf('%s.m',daesys);
    fid = fopen(daefile,'wt');
    fprintf(fid,'%s\n',strjoin({'function [res,ires] = ',daesys,'(t,y,ydot)'},''),'global UX',daeeqs{1:length(daeeqs)},'res = res(:);','ires = 0;','end');
    fclose(fid);

    for k=1:eqps.NumOfEqps
        for j=1:eqps.Equipment{k}.NumOfInitCond
            AddInitCondition(init,eqps.Equipment{k}.InitVars{j},eqps.Equipment{k}.InitValues(j))
        end
    end

    y0_fixed = zeros(size(y0,1),1); ydot0 = y0_fixed; ydot0_fixed = y0_fixed;
    for k=1:init.NumOfInitCond
        if strcmp(init.InitVars{k}(1:4),'diff')
            idx = GetVarIndex(vars,init.InitVars{k}(6:length(init.InitVars{k})-1));
            ydot0(vars.Variables{idx}.ID) = init.InitValues(k);
			ydot0_fixed(vars.Variables{idx}.ID) = true;
        else
            idx = GetVarIndex(vars,init.InitVars{k});
            if strcmp(vars.Variables{idx}.Name,init.InitVars{k})
                y0(vars.Variables{idx}.ID) = init.InitValues(k);
                y0_fixed(vars.Variables{idx}.ID) = true;
            end
        end
    end
end

function expr = ReplaceVariable(expr,vars)
    sMath = {' ',',','=','+','-','*','/','^','(',')','[',']','{','}'};
    for k=1:vars.NumOfVars
        varName = vars.Variables{k}.Name;
        sPos = strfind(expr,varName);
        varPos = [];
        for j=1:length(sPos)
            torf = true;
            if sPos(j) == 1  && torf
                for m=1:length(sMath)
                    if expr(sPos(j)+length(varName)) == sMath{m} && torf
                        varPos = [varPos,sPos(j)];
                        torf = false;
                    end
                end
            else
                if sPos(j)-1+length(varName) == length(expr) && torf
                    varPos = [varPos,sPos(j)];
                    torf = false;
                else
                    for m=1:length(sMath)
                        if expr(sPos(j)-1) == sMath{m} && torf
                            for n=1:length(sMath)
                                if expr(sPos(j)+length(varName)) == sMath{n} && torf
                                    varPos = [varPos,sPos(j)];
                                    torf = false;
                                end
                            end
                        end
                    end
                end      
            end
        end
        
        s = strcat('y(',num2str(vars.Variables{k}.ID),')');
        if length(varPos) ~= false
            if length(varPos) == 1
                if varPos == 1
                    expr = strjoin({s,expr(varPos+length(varName):length(expr))},'');
                elseif varPos == length(expr)
                    expr = strjoin({expr(1:varPos-1),s},'');
                else
                    expr = strjoin({expr(1:varPos-1),s,expr(varPos+length(varName):length(expr))},'');
                end
            else
                exprsplit = cell(length(varPos)+1,1);
                if varPos(1) == 1
                    exprsplit{1} = s;
                else
                    exprsplit{1} = strjoin({expr(1:varPos(1)-1),s},'');
                end
                if varPos(length(varPos)) == length(expr)
                    exprsplit{length(exprsplit)} = strjoin({expr(varPos(length(varPos))+length(varName):length(expr)),s},'');
                else
                    exprsplit{length(exprsplit)} = expr(varPos(length(varPos))+length(varName):length(expr));
                end
                for j=2:length(exprsplit)-1
                    exprsplit{j} = strjoin({expr(varPos(j-1)+length(varName):varPos(j)-1),s},'');
                end
                expr = strjoin(exprsplit,'');
            end
        end

        diffvar = strcat('diff(',s,')');
        sdiffPos = strfind(expr,diffvar);
        vardiffPos = [];
        for j=1:length(sdiffPos)
            torf = true;
            if sdiffPos(j) == 1  && torf
                for m=1:length(sMath)
                    if expr(sdiffPos(j)+length(diffvar)) == sMath{m} && torf
                        vardiffPos = [vardiffPos,sdiffPos(j)];
                        torf = false;
                    end
                end
            else
                if sdiffPos(j)-1+length(diffvar) == length(expr) && torf
                    vardiffPos = [vardiffPos,sdiffPos(j)];
                    torf = false;
                else
                    for m=1:length(sMath)
                        if expr(sdiffPos(j)-1) == sMath{m} && torf
                            for n=1:length(sMath)
                                if expr(sdiffPos(j)+length(diffvar)) == sMath{n} && torf
                                    vardiffPos = [vardiffPos,sdiffPos(j)];
                                    torf = false;
                                end
                            end
                        end
                    end
                end      
            end
        end
        
        sdiff = strcat('ydot(',num2str(vars.Variables{k}.ID),')');
        if length(vardiffPos) ~= false
            if length(vardiffPos) == 1
                if vardiffPos == 1
                    expr = strjoin({sdiff,expr(vardiffPos+length(diffvar):length(expr))},'');
                elseif vardiffPos == length(expr)
                    expr = strjoin({expr(1:vardiffPos-1),sdiff},'');
                else
                    expr = strjoin({expr(1:vardiffPos-1),sdiff,expr(vardiffPos+length(diffvar):length(expr))},'');
                end
            else
                exprsplit = cell(length(vardiffPos)+1,1);
                if vardiffPos(1) == 1
                    exprsplit{1} = sdiff;
                else
                    exprsplit{1} = strjoin({expr(1:vardiffPos(1)-1),sdiff},'');
                end
                if vardiffPos(length(vardiffPos)) == length(expr)
                    exprsplit{length(exprsplit)} = strjoin({expr(vardiffPos(length(vardiffPos))+length(diffvar):length(expr)),sdiff},'');
                else
                    exprsplit{length(exprsplit)} = expr(vardiffPos(length(vardiffPos))+length(diffvar):length(expr));
                end
                for j=2:length(exprsplit)-1
                    exprsplit{j} = strjoin({expr(vardiffPos(j-1)+length(diffvar):vardiffPos(j)-1),sdiff},'');
                end
                expr = strjoin(exprsplit,'');
            end
        end
    end
end
