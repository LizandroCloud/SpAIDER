function [y0,ydot0,resnorm] = findcic(resfile,y0,fixed_y0,ydot0,fixed_ydot0,atol)
%function [y0,ydot0,resnorm] = findcic(resfile,init,vars,y0_guess,ydot0_guess,atol)
%{
    initres = cell(init.NumOfInitCond,1);
    n = length(y0_guess);
    dbdy = zeros(init.NumOfInitCond,n);
    dbdydot = zeros(init.NumOfInitCond,n);
    for k=1:init.NumOfInitCond
        if strcmp(init.InitVars{k}(1:4),'diff')
            idx = GetVarIndex(vars,init.InitVars{k}(6:length(init.InitVars{k})-1));
            initres{k} = eval(strjoin({'@(y,ydot)ydot(',num2str(vars.Variables{idx}.ID),') - ',num2str(init.InitValues(k))},''));
            dbdydot(k,vars.Variables{idx}.ID) = 1;
        else
            idx = GetVarIndex(vars,init.InitVars{k});
            if strcmp(vars.Variables{idx}.Name,init.InitVars{k})
                initres{k} = eval(strjoin({'@(y,ydot)y(',num2str(vars.Variables{idx}.ID),') - ',num2str(init.InitValues(k))},''));
                dbdy(k,vars.Variables{idx}.ID) = 1;
            end
        end
    end
%}
    n = length(y0);
    if isempty(fixed_y0)
      free_y = 1:n;
      notfree_y = [];
    else
      free_y = find(fixed_y0 == 0);
      notfree_y = find(fixed_y0 == 1);
    end
    if isempty(fixed_ydot0)
      free_ydot = 1:n;
      notfree_ydot = [];
    else
      free_ydot = find(fixed_ydot0 == 0);
      notfree_ydot = find(fixed_ydot0 == 1);
    end
    if length(free_y) + length(free_ydot) < n
        error('Too many initial conditions specified');
    end
    
    % Eliminate from fmincon variables that are fixed
    ycic = [y0(free_y);ydot0(free_ydot)];
    
    % Send to residuals functions variables that are fixed.
    y_init = y0(notfree_y);
    ydot_init = ydot0(notfree_ydot);

    % Send to residuals functions fixed variables structure
    Cfree{1} = free_y;
    Cfree{2} = notfree_y;
    Cfree{3} = free_ydot;
    Cfree{4} = notfree_ydot;
    
    options = optimset(...
        'Algorithm','interior-point',...
        'Display','iter', ...
        'MaxIter',300, ...
        'MaxFunEvals',inf, ...
        'AlwaysHonorConstraints','none', ...
        'SubproblemAlgorithm','ldl-factorization',...
        'TolCon',atol, ...
        'TolFun',atol, ...
        'TolX',eps);
%    [ycic,~,exitflag] = fmincon(@(ycic)ResidualsNorm(resfile,ycic,initres),[y0_guess;ydot0_guess],[],[],[],[],[],[],@(ycic)ResidualsConstraint(resfile,ycic,initres),options);
%    y0 = ycic(1:length(ycic)/2);
%    ydot0 = ycic(length(ycic)/2+1:end);
    [ycic,resnorm] = fmincon(@(ycic)ResidualsNorm(resfile,ycic,y_init,ydot_init,Cfree),ycic,[],[],[],[],[],[],@(ycic)ResidualsConstraint(resfile,ycic,y_init,ydot_init,Cfree),options);
    y0(Cfree{1}) = ycic(1:length(Cfree{1}));
    y0(Cfree{2}) = y_init;
    ydot0(Cfree{3}) = ycic(length(Cfree{1})+1:end);
    ydot0(Cfree{4}) = ydot_init;


%{
    if exitflag == -1 || exitflag == -2  || exitflag == 0
        resf = daesys(0,y0,ydot0);
        resb = zeros(init.NumOfInitCond,1);
        for k=1:init.NumOfInitCond
            resb(k) = initres{k}(y0,ydot0);
        end
        resnorm = norm([resf;resb]);
        numOfiter = 300; k = 0;
        while resnorm > atol && k < numOfiter
            %Partial derivatives update for jacobian matrix
            daefun = @(t,y,ydot)resfile(t,y,ydot);
            [dfdy,~] = numjac(daefun,0,y0,resf,atol*ones(n,1),[],'off',[],[],ydot0);
            daefun = @(t,ydot,y)resfile(t,y,ydot);
            [dfdydot,~] = numjac(daefun,0,ydot0,resf,atol*ones(n,1),[],'off',[],[],y0);
            Jac = [[dfdy,dfdydot];[dbdy,dbdydot]];

            %QR decomposition to solve rectangular linear system
            %   0 = df + dfdy*dy + dfdydot*dydot
            %   0 = df + [dfdy,dfdydot]*[dy;dydot]

            %   0 = db + dbdy*dy + dbdydot*dydot
            %   0 = db + [dbdy,dbdydot]*[dy;dydot]

            %   0 = [df;db] + [[dfdy,dfydot];[dbdy,dbdydot]]*[dy;dydot]
            %   Jac = [[dfdy,dfydot];[dbdy,dbdydot]];
            %   [Q,R] = qr(Jac) --> [[dfdy,dfydot];[dbdy,dbdydot]] = Q*R
            %   0 = [df;db] + Q*R*[dy;dydot]
            %   0 = Q'*[df;db] + R*[dy;dydot] --> [dy;dydot] = R \ - Q'*[df;db]
            %   0 = [df;db] + Q*R*[dy;dydot]
            [Q,R,E] = qr(full(Jac)); r = rank(R);
            d = E * (R(1:r,:) \ (-Q(:,1:r)'*[resf;resb]));
            dy = d(1:n); dydot = d(n+1:end);
            y0 = y0 + dy; ydot0 = ydot0 + dydot;
            resf = daesys(0,y0,ydot0);
            resb = zeros(init.NumOfInitCond,1);
            for j=1:init.NumOfInitCond
                resb(j) = initres{j}(y0,ydot0);
            end
            resnorm = norm([resf;resb]);
            k=k+1;
        end
    end
%}
%{        
    y0 = y0_guess; ydot0 = ydot0_guess;
    resf = daesys(0,y0,ydot0,PP);
    resb = zeros(init.NumOfInitCond,1);
    for k=1:init.NumOfInitCond
        resb(k) = initres{k}(y0,ydot0);
    end
    resnorm = norm([resf;resb]);
    numOfiter = 300; k = 0;
    while resnorm > atol && k < numOfiter
        %Partial derivatives update for jacobian matrix
        [dfdy,~] = numjac(@(t,y,ydot,PP)daesys(t,y,ydot,PP),0,y0,resf,atol*ones(n,1),[],'off',[],[],ydot0,PP);
        [dfdydot,~] = numjac(@(t,ydot,y,PP)daesys(t,y,ydot,PP),0,ydot0,resf,atol*ones(n,1),[],'off',[],[],y0,PP);
        Jac = [[dfdy,dfdydot];[dbdy,dbdydot]];
        
        %QR decomposition to solve rectangular linear system
        %   0 = df + dfdy*dy + dfdydot*dydot
        %   0 = df + [dfdy,dfdydot]*[dy;dydot]

        %   0 = db + dbdy*dy + dbdydot*dydot
        %   0 = db + [dbdy,dbdydot]*[dy;dydot]

        %   0 = [df;db] + [[dfdy,dfydot];[dbdy,dbdydot]]*[dy;dydot]
        %   Jac = [[dfdy,dfydot];[dbdy,dbdydot]];
        %   [Q,R] = qr(Jac) --> [[dfdy,dfydot];[dbdy,dbdydot]] = Q*R
        %   0 = [df;db] + Q*R*[dy;dydot]
        %   0 = Q'*[df;db] + R*[dy;dydot] --> [dy;dydot] = R \ - Q'*[df;db]
        %   0 = [df;db] + Q*R*[dy;dydot]
        [Q,R,E] = qr(full(Jac)); r = rank(R);
        d = E * (R(1:r,:) \ (-Q(:,1:r)'*[resf;resb]));
        dy = d(1:n); dydot = d(n+1:end);
        y0 = y0 + dy; ydot0 = ydot0 + dydot;
        resf = daesys(0,y0,ydot0,PP);
        resb = zeros(init.NumOfInitCond,1);
        for j=1:init.NumOfInitCond
            resb(j) = initres{j}(y0,ydot0);
        end
        resnorm = norm([resf;resb]);
        k=k+1;
    end
%}
end

%{
function [y0,ydot0,exitflag] = findcic(init,vars,y0_guess,ydot0_guess,PP,atol)
    initres = cell(init.NumOfInitCond,1);
    for k=1:init.NumOfInitCond
        if strcmp(init.InitVars{k}(1:4),'diff')
            idx = GetVarIndex(vars,init.InitVars{k}(6:length(init.InitVars{k})-1));
            initres{k} = eval(strjoin({'@(y,ydot)ydot(',num2str(vars.Variables{idx}.ID),') - ',num2str(init.InitValues(k))},''));
        else
            idx = GetVarIndex(vars,init.InitVars{k});
            if strcmp(vars.Variables{idx}.Name,init.InitVars{k})
                initres{k} = eval(strjoin({'@(y,ydot)y(',num2str(vars.Variables{idx}.ID),') - ',num2str(init.InitValues(k))},''));
            end
        end
    end   
    options = optimset(...
        'Algorithm','interior-point',...
        'Display','iter', ...
        'MaxIter',300, ...
        'MaxFunEvals',inf, ...
        'AlwaysHonorConstraints','none', ...
        'SubproblemAlgorithm','ldl-factorization',...
        'TolCon',atol, ...
        'TolFun',atol, ...
        'TolX',eps);
    [ycic,~,exitflag] = fmincon(@(ycic)ResidualsNorm(ycic,initres,PP),[y0_guess;ydot0_guess],[],[],[],[],[],[],@(ycic)ResidualsConstaint(ycic,initres,PP),options);
    y0 = ycic(1:length(ycic)/2);
    ydot0 = ycic(length(ycic)/2+1:end);
%    yddot0_guess = zeros(length(ydot0_guess),1);
%    [ycic,resnrm] = fmincon(@(ycic)ResidualsNorm(ycic,initres,PP),[y0_guess;ydot0_guess;yddot0_guess],[],[],[],[],[],[],@(ycic)ResidualsConstaint(ycic,initres,PP),options);
%    y0 = ycic(1:length(ycic)/3);
%    ydot0 = ycic(length(ycic)/3+1:length(ycic)*(2/3));
end
%}

%function f = ResidualsNorm(resfile,ycic,initres)
function f = ResidualsNorm(resfile,ycic,y_init,ydot_init,Cfree)
%    [~,res] = ResidualsConstraint(resfile,ycic,initres);
    [~,res] = ResidualsConstraint(resfile,ycic,y_init,ydot_init,Cfree);
    f = norm(res);
end

%function [resineq,reseq] = ResidualsConstraint(resfile,ycic,initres)
function [resineq,reseq] = ResidualsConstraint(resfile,ycic,y_init,ydot_init,Cfree)
    y = zeros(length(Cfree{1}) + length(Cfree{2}),1);
    y(Cfree{1}) = ycic(1:length(Cfree{1}));
    y(Cfree{2}) = y_init;

    ydot = zeros(length(Cfree{3}) + length(Cfree{4}),1);
    ydot(Cfree{3}) = ycic(length(Cfree{1})+1:end);
    ydot(Cfree{4}) = ydot_init;
    
    reseq = feval(resfile,0,y,ydot);
    
%    n = length(ycic)/2;
%    reseq = feval(resfile,0,ycic(1:n),ycic(n+1:end));
%    for k=1:length(initres)
%        reseq(n+k) = initres{k}(ycic(1:n),ycic(n+1:end));
%    end	
%{
    n = length(ycic)/3;
    reseq = daesys(0,ycic(1:n),ycic(n+1:2*n),PP);
    h = 1e-8;
    dreseqdt = (1/h)*( ...
        daesys(h,ycic(1:n) + h*ycic(n+1:2*n),ycic(n+1:2*n) + h*ycic(2*n+1:end),PP) - reseq);
	reseq = [reseq;dreseqdt];
    for k=1:length(initres)
        reseq(2*n+k) = initres{k}(ycic(1:n),ycic(n+1:2*n));
    end
%}
    reseq=reseq(:);
    resineq = []; 
end
