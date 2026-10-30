classdef PID_pos < ModelTypes
	properties
		action = {};	%Controller action. Valid values = {'on','off'}
	end
	
	methods
        function obj = PID_pos(Name,varargin)
        %MODEL BRIEF
            obj.Description = 'Proportional-Integral-Derivative control positional law';
            obj.Name = Name;
			
        %PARAMETERS
        %User-defined parameters
            if length(varargin) > 2
                error('Too many input parameters as arguments. Valid parameters: action');
            elseif length(varargin) > 1
                k=1;
                while k<length(varargin)
                    if strcmp(varargin{k},'action')
                        obj.action = varargin{k+1};
                    k=k+1;
                    end
                end
            end

        %VARIABLES
            AddVariable(obj,'PV',admensional);
            AddVariable(obj,'SP',admensional);
            AddVariable(obj,'MVoutput',admensional);
            AddVariable(obj,'MVinput',admensional);
            AddVariable(obj,'unboutp',admensional);
            AddVariable(obj,'errorC',admensional);
            AddVariable(obj,'derivFilt',admensional);
            AddVariable(obj,'AWBTintterm',admensional);
            AddVariable(obj,'KC',admensional);
            AddVariable(obj,'TI',time_sec);
			AddVariable(obj,'KI',frequency);
            AddVariable(obj,'TD',time_sec);
			AddVariable(obj,'KD',time_sec);
			AddVariable(obj,'mode',admensional);
			

            idx = GetVarIndex(obj,'PV');    		obj.Variables{idx}.Name = strcat(Name,'.PV');   		obj.Variables{idx}.Description = 'Bounded input signal from transmitter';   				obj.Variables{idx}.Value = 0.5;		obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'SP');    		obj.Variables{idx}.Name = strcat(Name,'.SP');   		obj.Variables{idx}.Description = 'Bounded set-point from transmitter';      				obj.Variables{idx}.Value = 0.5; 	obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'MVoutput');  	obj.Variables{idx}.Name = strcat(Name,'.MVoutput');   	obj.Variables{idx}.Description = 'Bounded manipulated variable output signal';   			obj.Variables{idx}.Value = 0.5; 	obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'MVinput');   	obj.Variables{idx}.Name = strcat(Name,'.MVinput');   	obj.Variables{idx}.Description = 'Manipulated variable input signal (manual mode)'; 		obj.Variables{idx}.Value = 0.5; 	obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'unboutp');       obj.Variables{idx}.Name = strcat(Name,'.unboutp');      obj.Variables{idx}.Description = 'Unbounded output signal';           						obj.Variables{idx}.Value = 0.5; 	obj.Variables{idx}.Unit = '';
			idx = GetVarIndex(obj,'errorC');       	obj.Variables{idx}.Name = strcat(Name,'.errorC');      	obj.Variables{idx}.Description = 'Controller error';           								obj.Variables{idx}.Value = 0; 		obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'derivFilt');  	obj.Variables{idx}.Name = strcat(Name,'.derivFilt'); 	obj.Variables{idx}.Description = 'Derivative term filter';                                  obj.Variables{idx}.Value = 0;  		obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'AWBTintterm');  	obj.Variables{idx}.Name = strcat(Name,'.AWBTintterm'); 	obj.Variables{idx}.Description = 'Integral term with anti-windup and bumpless transfer';   	obj.Variables{idx}.Value = 0;  		obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'KC');     		obj.Variables{idx}.Name = strcat(Name,'.KC');    		obj.Variables{idx}.Description = 'Proportional gain';             							obj.Variables{idx}.Value = 1; 		obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'TI');         	obj.Variables{idx}.Name = strcat(Name,'.TI');        	obj.Variables{idx}.Description = 'Integral time';                          					obj.Variables{idx}.Value = 10; 		obj.Variables{idx}.Unit = 's';
            idx = GetVarIndex(obj,'KI');     		obj.Variables{idx}.Name = strcat(Name,'.KI');    		obj.Variables{idx}.Description = 'Integral gain';                            				obj.Variables{idx}.Value = 1/10;	obj.Variables{idx}.Unit = 'Hz';
			idx = GetVarIndex(obj,'TD');         	obj.Variables{idx}.Name = strcat(Name,'.TD');        	obj.Variables{idx}.Description = 'Derivative time';                          				obj.Variables{idx}.Value = 0; 		obj.Variables{idx}.Unit = 's';
            idx = GetVarIndex(obj,'KD');     		obj.Variables{idx}.Name = strcat(Name,'.KD');    		obj.Variables{idx}.Description = 'Derivative gain';                            				obj.Variables{idx}.Value = 0; 		obj.Variables{idx}.Unit = 's';
            idx = GetVarIndex(obj,'mode');     		obj.Variables{idx}.Name = strcat(Name,'.mode');    		obj.Variables{idx}.Description = 'Controler mode: Manual=1, Auto=0';                        obj.Variables{idx}.Value = 0; 		obj.Variables{idx}.Unit = '';			

        %Variable tokens
            PV			= strcat(obj.Name,'.PV');
            SP          = strcat(obj.Name,'.SP');
            MVoutput    = strcat(obj.Name,'.MVoutput');
            MVinput     = strcat(obj.Name,'.MVinput');
            unboutp     = strcat(obj.Name,'.unboutp');
            errorC      = strcat(obj.Name,'.errorC');
            derivFilt      = strcat(obj.Name,'.derivFilt');
            AWBTintterm	= strcat(obj.Name,'.AWBTintterm');
            KC          = strcat(obj.Name,'.KC');
            TI          = strcat(obj.Name,'.TI');
			KI          = strcat(obj.Name,'.KI');
			TD          = strcat(obj.Name,'.TD');
			KD          = strcat(obj.Name,'.KD');
			mode        = strcat(obj.Name,'.mode');

        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);

        %EQUATIONS
        %Controller error
			AddEquation(obj,strcat(errorC,' -( ',SP,' - ',PV,')'));

        %Derivative term filter
			AddEquation(obj,strcat('0.1*(',TD,'+1e-3)*diff(',derivFilt,')',' -( ',errorC,' - ',derivFilt,')'));

        %Integral term with anti-windup and bumpless transfer
			AddEquation(obj,strcat('diff(',AWBTintterm,')',' -( ',KI,'*',errorC,'*', ...
				'0.5*(1+tanh((',unboutp,'-1e-5)/1e-6))*0.5*(1+tanh(((1-1e-5)-',unboutp,')/1e-6))',')'));

        %Conventional x decoupled PID law: integral action parameter
			AddEquation(obj,strcat(KI,'*',TI,' - ',KC));

        %Conventional x decoupled PID law: derivative action parameter
			AddEquation(obj,strcat(KD,' - ',KC,'*',TD));

        %Bounded output
			AddEquation(obj,strcat(MVoutput,' -( ','0.5*(1+tanh((2*',mode,'-1)/1e-5))*',MVinput, ...
				' + 0.5*(1+tanh((1-2*',mode,')/1e-5))*(', ...
					'0.5*(1+tanh((',unboutp,'-1e-5)/1e-4))*0.5*(1+tanh(((1-1e-5)-',unboutp,')/1e-4))*',unboutp,...
					' + 0.5*(1+tanh((',unboutp,'-(1-1e-5))/1e-4))',')', ...
				')'));

            switch obj.action
                case 'Direct'
                %PID direct positional law
                    AddEquation(obj,strcat(unboutp,' + ','(',KC,'*',errorC,' + ',AWBTintterm,' + ',KD,'*diff(',derivFilt,')',')'));

                case 'Reverse'
                %PID reverse positional law
                    AddEquation(obj,strcat(unboutp,' - ','(',KC,'*',errorC,' + ',AWBTintterm,' + ',KD,'*diff(',derivFilt,')',')'));

                otherwise
                    error('invalid controller action mode. Valid inputs: {Direct,Reverse}');
            end

        %INITIAL CONDITIONS
			AddInitCondition(obj,strcat('diff(',derivFilt,')'),0.0);
        end
	end
end