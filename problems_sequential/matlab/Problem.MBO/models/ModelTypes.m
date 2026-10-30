classdef ModelTypes < handle
    properties
        Name = ''
        Description = '';
        Parameters = {};
        NumOfPars = 0;
        Variables = {};
        NumOfVars = 0;
        Equipment = {};
        NumOfEqps = 0;
        Equations = {};
        NumOfEqs = 0;
        InitVars = {};
        InitValues = [];
        NumOfInitCond = 0;
    end
    
    methods
 %{
        function [] = AddParameter(obj,ParName,ParType)
            SetParameter(obj,ParName,ParType);
        end
        
        function [] = SetParameter(obj,ParName,ParType)
            k = GetNumOfPars(obj);
            SetNumOfPars(obj,k+1);
            obj.Parameters{k+1,1} = ParType;
            obj.Parameters{k+1}.Name = ParName;
        end
        
        function obj = SetNumOfPars(obj,value)
           obj.NumOfPars = value;
        end
        
        function value = GetNumOfPars(obj)
           value = obj.NumOfPars;
        end
        
        function idx = GetParIndex(obj,ParName)
            torf = true;
            for k=1:length(obj.Parameters)
                if strcmp(obj.Parameters{k}.Name,ParName) == true && torf
                    idx = k;
                    torf = false;
                end
            end
        end

        function [] = SetParameterValue(obj,ParName,Value)
            idx = GetParIndex(obj,ParName);
            obj.Parameters{idx}.Value = Value;
        end
        
        function Value = GetParameterValue(obj,ParName)
            idx = GetParIndex(obj,ParName);
            Value = obj.Parameters{idx}.Value;
        end
        
        function obj = SystemParameters(obj)
            for k=1:obj.NumOfEqps
                for j=1:obj.Equipment{k}.NumOfPars
                    AddParameter(obj,obj.Equipment{k}.Parameters{j}.Name,obj.Equipment{k}.Parameters{j});
                end
            end
        end
%}        
        function [] = AddVariable(obj,VarName,VarType)
            SetVariable(obj,VarName,VarType);
        end
        
        function [] = SetVariable(obj,VarName,VarType)
            k = GetNumOfVars(obj);
            SetNumOfVars(obj,k+1);
            obj.Variables{k+1,1} = VarType;
            obj.Variables{k+1}.Name = VarName;
        end
        
        function obj = SetNumOfVars(obj,value)
           obj.NumOfVars = value;
        end
        
        function value = GetNumOfVars(obj)
           value = obj.NumOfVars;
        end
        
        function idx = GetVarIndex(obj,VarName)
            torf = true;
            for k=1:length(obj.Variables)
                if strcmp(obj.Variables{k}.Name,VarName) == true && torf
                    idx = k;
                    torf = false;
                end
            end
        end

        function obj = SystemVariables(obj)
            for k=1:obj.NumOfEqps
                for j=1:obj.Equipment{k}.NumOfVars
                    AddVariable(obj,obj.Equipment{k}.Variables{j}.Name,obj.Equipment{k}.Variables{j});
                end
            end
        end
        
        function [] = AddEquipment(obj,EqpType)
            SetEquipment(obj,EqpType);
        end
        
        function [] = SetEquipment(obj,EqpType)
            k = GetNumOfEqps(obj);
            SetNumOfEqps(obj,k+1);
            obj.Equipment{k+1,1} = EqpType;
        end
        
        function obj = SetNumOfEqps(obj,value)
           obj.NumOfEqps = value;
        end
        
        function value = GetNumOfEqps(obj)
           value = obj.NumOfEqps;
        end

        function idx = GetEqpIndex(obj,EqpName)
            torf = true;
            for k=1:obj.NumOfEqps
                if strcmp(obj.Equipment{k}.Name,EqpName) == true && torf
                    idx = k;
                    torf = false;
                end
            end
        end
        
        function [] = AddEquation(obj,expr)
            SetEquation(obj,expr);
        end
        
        function [] = SetEquation(obj,expr)
            k = GetNumOfEqs(obj);
            SetNumOfEqs(obj,k+1);
            obj.Equations{k+1,1} = expr;
        end
        
        function obj = SetNumOfEqs(obj,value)
           obj.NumOfEqs = value;
        end
        
        function value = GetNumOfEqs(obj)
           value = obj.NumOfEqs;
        end
        
        function obj = SystemEquations(obj)
            for k=1:obj.NumOfEqps
                for j=1:obj.Equipment{k}.NumOfEqs
                    AddEquation(obj,obj.Equipment{k}.Equations{j});
                end
            end
        end
        
        function [] = AddInitCondition(obj,VarName,InitValue)
            SetInitCondition(obj,VarName,InitValue);
        end
        
        function [] = SetInitCondition(obj,VarName,InitValue)
            k = GetNumOfInitCond(obj);
            SetNumOfInitCond(obj,k+1);
            obj.InitValues(k+1,1) = InitValue;
            obj.InitVars{k+1,1} = VarName;
        end
        
        function obj = SetNumOfInitCond(obj,value)
           obj.NumOfInitCond = value;
        end
        
        function value = GetNumOfInitCond(obj)
           value = obj.NumOfInitCond;
        end
        
        function obj = SystemInitConditions(obj)
            for k=1:obj.NumOfEqps
                for j=1:obj.Equipment{k}.NumOfInitCond
                    AddInitCondition(obj,obj.Equipment{k}.InitVars{j},obj.Equipment{k}.InitValues(j))
                end
            end
        end
    end
end

