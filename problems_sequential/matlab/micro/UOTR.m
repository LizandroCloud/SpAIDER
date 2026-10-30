function [pfd, y0,ydot0] = UOTR(~)
    %PHYSICAL PROPERTIES
    global PP
    PP = mbo;

    %VARIABLES
    vars = variables_head;
    
    %EQUIPMENT
    eqps = equipment_head;
    AddEquipment(eqps,triphase_source('feed'));
 %   AddEquipment(eqps,triphase_sink('liq'));
    AddEquipment(eqps,gasliq_sink('gas'));
	AddEquipment(eqps,water_sink('water'));
    
    AddEquipment(eqps,heater_liq('P_002'));
    AddEquipment(eqps,degasser('V_TO_001','L',8.00,'D',2.50));
	AddEquipment(eqps,oil_treater('TO_001','L',15.0,'D',4.0,'TRL50',210,'b',3.0));

    AddEquipment(eqps,temperature_transmitter('TIT_001','RangeMax',(120+273.15),'RangeMin',(0+273.15)));
    AddEquipment(eqps,PID_pos('TIC_001','action','Reverse'));    
    
    AddEquipment(eqps,pressure_transmitter('PIT_002','RangeMax',1500,'RangeMin',0.00));
    AddEquipment(eqps,PID_pos('PIC_002','action','Direct'));

%    AddEquipment(eqps,valv_liq_ctrl('LV_003','Cv',1800,'FL',0.91));
    AddEquipment(eqps,level_transmitter('LIT_003','RangeMax',2.50,'RangeMin',0.00));
    AddEquipment(eqps,PID_pos('LIC_003','action','Direct'));
 
    AddEquipment(eqps,level_transmitter('LIT_004','RangeMax',4.00,'RangeMin',0.00));
    AddEquipment(eqps,PID_pos('LIC_004','action','Direct'));
 
    %CONNECTIONS
    cnct = connections_head;
    AddEqpConnection(cnct,'feed.Outlet','P_002.Inlet');
%    AddEqpConnection(cnct,'LV_003.Outlet','liq.Inlet');
%    AddEqpConnection(cnct,'TO_001.OutletL','liq.Inlet');
    AddEqpConnection(cnct,'V_TO_001.OutletG','gas.Inlet');
    AddEqpConnection(cnct,'TO_001.OutletW','water.Inlet');
    
    AddEqpConnection(cnct,'P_002.Outlet','V_TO_001.Inlet');
	AddEqpConnection(cnct,'V_TO_001.OutletL','TO_001.Inlet');
%    AddEqpConnection(cnct,'V_TO_001.OutletG','PV_002.Inlet');
    
%	 AddEqpConnection(cnct,'TO_001.OutletL','LV_003.Inlet');
%    AddEqpConnection(cnct,'TO_001.OutletW','LV_004.Inlet');

    AddVarConnection(cnct,'P_002.Outlet.T','TIT_001.PVinput');
    AddVarConnection(cnct,'TIT_001.PVoutput','TIC_001.PV');
    AddVarConnection(cnct,'TIT_001.SPoutput','TIC_001.SP');    
    
    AddVarConnection(cnct,'V_TO_001.OutletG.P','PIT_002.PVinput');
    AddVarConnection(cnct,'PIT_002.PVoutput','PIC_002.PV');
    AddVarConnection(cnct,'PIT_002.SPoutput','PIC_002.SP');
    
    AddVarConnection(cnct,'V_TO_001.Level','LIT_003.PVinput');
    AddVarConnection(cnct,'LIT_003.PVoutput','LIC_003.PV');
    AddVarConnection(cnct,'LIT_003.SPoutput','LIC_003.SP');
 %   AddVarConnection(cnct,'LIC_003.MVoutput','LV_003.MVinput');

    AddVarConnection(cnct,'TO_001.LevelW','LIT_004.PVinput');
    AddVarConnection(cnct,'LIT_004.PVoutput','LIC_004.PV');
    AddVarConnection(cnct,'LIT_004.SPoutput','LIC_004.SP');
    
    %SPECIFICATIONS
    spec = specifications_head;
    AddSpecification(spec,'feed.Outlet.Qf - 22000');
    AddSpecification(spec,'feed.Outlet.Qg - feed.Outlet.Qf*5.00');
    AddSpecification(spec,'feed.Outlet.BSW - 0.10');
    AddSpecification(spec,'feed.Outlet.T - (65+273.15)');
    AddSpecification(spec,'feed.Outlet.d60 - 141.5/(30+131.5)');
    AddSpecification(spec,'feed.Outlet.Tbn - (11.8*feed.Outlet.d60)^3/1.8');
    AddSpecification(spec,'feed.Outlet.dg - 1.05');
    AddSpecification(spec,'feed.Outlet.Ycg - 5.5');
    AddSpecification(spec,'feed.Outlet.ws - (255/PP.rho_H2Ostd)/(1 + (255/PP.rho_H2Ostd))');

    AddSpecification(spec,'TIT_001.SPinput - (80+273.15)');
    AddSpecification(spec,'TIC_001.KC - 1.00');
    AddSpecification(spec,'TIC_001.TI - 60');
    AddSpecification(spec,'TIC_001.TD - 0.00');
    AddSpecification(spec,'TIC_001.mode - 0.00');
    AddSpecification(spec,'TIC_001.MVinput - TIC_001.MVoutput');    
    
    AddSpecification(spec,'PIT_002.SPinput - (700+101.325)');
    AddSpecification(spec,'PIC_002.KC - 4.00');
    AddSpecification(spec,'PIC_002.TI - 60');
    AddSpecification(spec,'PIC_002.TD - 0.00');
    AddSpecification(spec,'PIC_002.mode - 0.00');
    AddSpecification(spec,'PIC_002.MVinput - PIC_002.MVoutput');

    AddSpecification(spec,'LIT_003.SPinput - 1.00');
    AddSpecification(spec,'LIC_003.KC - 5.00');
    AddSpecification(spec,'LIC_003.TI - 600');
    AddSpecification(spec,'LIC_003.TD - 0.00');
    AddSpecification(spec,'LIC_003.mode - 0.00');
    AddSpecification(spec,'LIC_003.MVinput - LIC_003.MVoutput');

    AddSpecification(spec,'LIT_004.SPinput - 0.5');
    AddSpecification(spec,'LIC_004.KC - 4.00');
    AddSpecification(spec,'LIC_004.TI - 600');
    AddSpecification(spec,'LIC_004.TD - 0.00');
    AddSpecification(spec,'LIC_004.mode - 0.00');
    AddSpecification(spec,'LIC_004.MVinput - LIC_004.MVoutput');
    
    %EQUATIONS
    eqs = equations_head;
    AddEquation(eqs,'P_002.DP - 120');
 %   AddEquation(eqs,'LV_003.Outlet.P - (200+101.325)');
    
    AddEquation(eqs,'P_002.Qin - (0 + TIC_001.MVoutput*(50e3 - 0))');
    AddEquation(eqs,'gas.Inlet.Qg - PIC_002.MVoutput*5e5');
%    AddEquation(eqs,'liq.Inlet.Qf - LIC_003.MVoutput*30000');
    AddEquation(eqs,'TO_001.OutletL.Qf - LIC_003.MVoutput*30000');
	AddEquation(eqs,'water.Inlet.Qw - LIC_004.MVoutput*20000');
    
    %INITIAL CONDITIONS
    init = initial_head;
    AddInitCondition(init,'V_TO_001.OutletG.P',(800+101.325));
    AddInitCondition(init,'V_TO_001.OutletG.T',(80+273.15));
    AddInitCondition(init,'V_TO_001.OutletL.BSW',0.20);
    AddInitCondition(init,'V_TO_001.OutletL.d60',141.5/(29+131.5));
    AddInitCondition(init,'V_TO_001.OutletL.Tbn',(11.3*141.5/(29+131.5))^3/1.8);
    AddInitCondition(init,'V_TO_001.OutletG.dg',0.7);
    AddInitCondition(init,'V_TO_001.OutletG.Ycg',4.0);
    AddInitCondition(init,'V_TO_001.OutletL.ws',(255/PP.rho_H2Ostd)/(1 + (255/PP.rho_H2Ostd)));
    AddInitCondition(init,'V_TO_001.Level',1.0);

    AddInitCondition(init,'TO_001.OutletW.T',(80+273.15));
    AddInitCondition(init,'TO_001.OutletL.d60',141.5/(29+131.5));
    AddInitCondition(init,'TO_001.OutletL.Tbn',(11.3*141.5/(29+131.5))^3/1.8);
    AddInitCondition(init,'TO_001.OutletW.ws',(255/PP.rho_H2Ostd)/(1 + (255/PP.rho_H2Ostd)));
    AddInitCondition(init,'TO_001.LevelW',1.0);

%    AddInitCondition(init,'P_002.Qin',10e3);
    
%    AddInitCondition(init,'LV_003.c',22/30);
%    AddInitCondition(init,'LIC_003.MVoutput',0.5);
    
    AddInitCondition(init,'TIC_001.MVoutput',25/50);
    AddInitCondition(init,'PIC_002.MVoutput',110/500);
    AddInitCondition(init,'LIC_003.MVoutput',22*(1-0.2)/30);
	AddInitCondition(init,'LIC_004.MVoutput',22*0.2/20);	
    
    %DAE SYSTEM BUILD
    tic
    [y0,y0_fixed,ydot0,ydot0_fixed] = BuildDAEsystem('daesys',vars,eqps,cnct,spec,eqs,init);
    elapsedTime = toc;
    disp(strjoin({'DAE system was built successfully in',num2str(elapsedTime),'seconds'}));

    %CONSISTENT INITIAL CONDITIONS SOLVE
    tic
    atol=1e-8;
    [y0,ydot0,resnorm] = findcic('daesys',y0,y0_fixed,ydot0,ydot0_fixed,atol);
    elapsedTime = toc;
    if resnorm < atol
        disp(strjoin({'Consistent initial conditions were obtained successfully in',num2str(elapsedTime),'seconds'}));
    else
       error(strjoin({'non-linear solver could not obtain consistent initial conditions. resnorm =',num2str(resnorm)}));
    end
%{    
    tic
    atol=1e-8;
    [y0,ydot0,exitflag] = findcic(init,vars,y0,ydot0,PP,atol);
    elapsedTime = toc;
    switch exitflag
        case -1
            error(strjoin({'non-linear solver could not obtain consistent initial conditions. An error on residuals function was found'}));
        case -2
            error(strjoin({'non-linear solver could not obtain consistent initial conditions. No feasible point was found.'}));
        case 0
            error(strjoin({'non-linear solver could not obtain consistent initial conditions. Number of iterations exceeded options.MaxIter'}));
        case 1
            disp(strjoin({'Consistent initial conditions were obtained successfully in',num2str(elapsedTime),'seconds'}));
        case 2
            disp(strjoin({'Consistent initial conditions were obtained successfully in',num2str(elapsedTime),'seconds'}));
            disp(strjoin({'Change in y0 was less than options.TolX and residuals function was less than',num2str(atol)}));
        case -3
            disp(strjoin({'Consistent initial conditions were obtained successfully in',num2str(elapsedTime),'seconds'}));
            disp(strjoin({'Residuals function at current iteration went below',num2str(atol)}));
    end
%}
%{
    %CONSISTENT INITIAL CONDITIONS SOLVE
    tic
    n = length(y0); atol=1e-8;
    cicopt = optimset('Algorithm','levenberg-marquardt', ...
                      'ScaleProblem','jacobian', ...
                      'InitDamping',eps, ...
                      'Display','iter', ...
                      'MaxFunEvals',inf, ...
                      'MaxIter',300, ...
                      'TolFun',atol, ...
                      'TolX',atol);
    [y0,res] = fsolve(@(y)daesys(0,y,zeros(n,1),PP),y0,cicopt);
    elapsedTime = toc;
    if norm(res) < atol
        disp(strjoin({'Steady-state conditions were obtained successfully in',num2str(elapsedTime),'seconds'}));
    else
        error(strjoin({'non-linear solver could not obtain steady-state conditions. resnorm =',num2str(norm(res))}));
    end
%}
	%DAE SYSTEM SOLVE
    tspan = (0:1:10800); rtol=1e-5; atol=1e-8;
    daeopt = odeset('RelTol',rtol, ...
                    'AbsTol',atol, ...
                    'NormControl','on', ...
                    'stats','on');
    tic
%    [pfd.time,y] = ode15i(@(t,y,ydot)daesys(t,y,ydot,PP),tspan,y0,ydot0,daeopt);
    [pfd.time,y] = ode15i(@(t,y,ydot)daesys(t,y,ydot),tspan,y0,ydot0,daeopt);
%    [pfd.time,y] = dasslc('daesys',tspan,y0,ydot0,PP,rtol,atol);
    elapsedTime = toc;
    disp(strjoin({'Simulation has finished successfully in',num2str(elapsedTime),'seconds'}));
    
    %DAE SYSTEM DECRYPT
    pfd = DecryptSystem(y,vars,eqps,pfd);
end

