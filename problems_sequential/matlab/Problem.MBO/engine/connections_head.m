classdef connections_head < handle
    properties
        EqpNameFrom = {};
        EqpNameTo = {};
        NumOfEqpCnct = 0;
        VarNameFrom = {};
        VarNameTo = {};
        NumOfVarCnct = 0;
    end
    
    methods
        function [] = AddEqpConnection(obj,EqpNameFrom,EqpNameTo)
            SetEqpConnection(obj,EqpNameFrom,EqpNameTo);
        end
        
        function [] = SetEqpConnection(obj,EqpNameFrom,EqpNameTo)
            k = GetNumOfEqpCnct(obj);
            SetNumOfEqpCnct(obj,k+1);
            obj.EqpNameFrom{k+1,1} = EqpNameFrom;
            obj.EqpNameTo{k+1,1} = EqpNameTo;
        end
        
        function obj = SetNumOfEqpCnct(obj,value)
           obj.NumOfEqpCnct = value;
        end
        
        function value = GetNumOfEqpCnct(obj)
           value = obj.NumOfEqpCnct;
        end
        
        function [] = AddVarConnection(obj,VarNameFrom,VarNameTo)
            SetVarConnection(obj,VarNameFrom,VarNameTo);
        end
        
        function [] = SetVarConnection(obj,VarNameFrom,VarNameTo)
            k = GetNumOfVarCnct(obj);
            SetNumOfVarCnct(obj,k+1);
            obj.VarNameFrom{k+1,1} = VarNameFrom;
            obj.VarNameTo{k+1,1} = VarNameTo;
        end
        
        function obj = SetNumOfVarCnct(obj,value)
           obj.NumOfVarCnct = value;
        end
        
        function value = GetNumOfVarCnct(obj)
           value = obj.NumOfVarCnct;
        end        
    end

end

