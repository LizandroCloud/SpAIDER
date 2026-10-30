classdef valv_liq_ctrl < valv_liq  
    properties
		tau	= 0;		%Valve time constant					[s]
	end
	
	methods
        function obj = valv_liq_ctrl(Name,varargin)
        %MODEL BRIEF
            obj = obj@valv_liq(Name,varargin{:});
			obj.Description = 'Basic equations for liquid stream flowing through valves for control applications';
            obj.Name = Name;
            
        %PARAMETERS
        %User-defined parameters
            if length(varargin) > 6
                error('Too many input parameters as arguments. Valid parameters: {Cv,FL,tau}');
            elseif length(varargin) > 1
                k=1;
                while k<length(varargin)
					switch varargin{k}
%						case 'Cv'
%							obj.Cv = varargin{k+1};
%						case 'FL'
%							obj.FL = varargin{k+1};
						case 'tau'
							obj.tau = varargin{k+1};
					end
                    k=k+1;
                end
            end
            
        %Dependent parameters       
%            obj.Kv = obj.Cv*(0.227124707)*sqrt(14.696/1.01325);
%			obj.Av = obj.Cv*(0.0000630902)*sqrt(999.054*14.696/101325);

        %Parameter tokens
			taupar = num2str(obj.tau);            
            
        %VARIABLES
            AddVariable(obj,'MVinput',admensional);

            idx = GetVarIndex(obj,'MVinput');	obj.Variables{idx}.Name = strcat(Name,'.MVinput');	obj.Variables{idx}.Description = 'Input signal from controller';	obj.Variables{idx}.Value = 0.5;	obj.Variables{idx}.Unit = '';

        %Variable tokens
            MVinput          	= strcat(obj.Name,'.MVinput');
			c          			= strcat(obj.Name,'.c');    
            
        %REVERSE MODEL INHERITANCE
%            SystemVariables(obj);
%            SystemEquations(obj);
%            SystemInitConditions(obj);
        
        %EQUATIONS
        %Linear first-order dynamic behavior
            AddEquation(obj,strcat('(1e-3 +',taupar,')*diff(',c,') + ',c,' - ',MVinput));

        end
    end
    
end

