classdef temperature_transmitter < ModelTypes
	properties
		RangeMax = temperature;	%Transmitter maximum range
		RangeMin = temperature;	%Transmitter minimum range
	end
	
	methods
        function obj = temperature_transmitter(Name,varargin)
        %MODEL BRIEF
            obj.Description = 'Transmitter model for temperature variables';
            obj.Name = Name;
			
        %PARAMETERS
        %User-defined parameters
            if length(varargin) > 4
                error('Too many input parameters as arguments. Valid parameters: {RangeMax,RangeMin}');
            elseif length(varargin) > 1
                k=1;
                while k<length(varargin)
                    if strcmp(varargin{k},'RangeMax')
                        obj.RangeMax.Value = varargin{k+1};
                    elseif strcmp(varargin{k},'RangeMin')
						obj.RangeMin.Value = varargin{k+1};
                    end
                    k=k+1;
                end
            end

        %Parameter tokens
            RangeMaxPar = num2str(obj.RangeMax.Value);
			RangeMinPar = num2str(obj.RangeMin.Value);

        %VARIABLES
            AddVariable(obj,'PVinput',temperature);
            AddVariable(obj,'SPinput',temperature);
            AddVariable(obj,'PVoutput',admensional);
            AddVariable(obj,'SPoutput',admensional);
			

            idx = GetVarIndex(obj,'PVinput');	obj.Variables{idx}.Name = strcat(Name,'.PVinput');	obj.Variables{idx}.Description = 'Input signal from process variable';   	obj.Variables{idx}.Value = 5000;	obj.Variables{idx}.Unit = 'K';
            idx = GetVarIndex(obj,'SPinput');   obj.Variables{idx}.Name = strcat(Name,'.SPinput');  obj.Variables{idx}.Description = 'Set-point from user/master loop';      	obj.Variables{idx}.Value = 5000; 	obj.Variables{idx}.Unit = 'K';
            idx = GetVarIndex(obj,'PVoutput');  obj.Variables{idx}.Name = strcat(Name,'.PVoutput'); obj.Variables{idx}.Description = 'Bounded process variable output signal';	obj.Variables{idx}.Value = 0.5; 	obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'SPoutput');  obj.Variables{idx}.Name = strcat(Name,'.SPoutput'); obj.Variables{idx}.Description = 'Bounded set-point output signal'; 		obj.Variables{idx}.Value = 0.5; 	obj.Variables{idx}.Unit = '';

        %Variable tokens
            PVinput		= strcat(obj.Name,'.PVinput');
            SPinput		= strcat(obj.Name,'.SPinput');
            PVoutput	= strcat(obj.Name,'.PVoutput');
            SPoutput    = strcat(obj.Name,'.SPoutput');

        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);

        %EQUATIONS
        %Bounded process variable output signal to controller input
			AddEquation(obj,strcat(PVoutput,'*(',RangeMaxPar,'-',RangeMinPar,')',' -( ',PVinput,' - ',RangeMinPar,')'));

        %Bounded set-point output signal to controller input
			AddEquation(obj,strcat(SPoutput,'*(',RangeMaxPar,'-',RangeMinPar,')',' -( ',SPinput,' - ',RangeMinPar,')'));

		
		end
	end
end

