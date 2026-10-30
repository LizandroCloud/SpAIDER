classdef blackoil
    properties (Constant)
        Pref            = 101.325               %Ref. pressure for std. cond.       [kPa]
        Tref            = (60-32)/1.8+273.15    %Ref. temperature for std. cond.    [K]
        R               = 8.31447               %Ideal gas universal constant       [kJ/kmol K]
        rho_H2Ostd      = 999.054               %Pure water mass density            [kg/m^3]
        rho_AIRstd      = 1.2938                %Air mass density                   [kg/m^3]
        MM_Air          = 29.0                  %Air molecular weight               [kg/kmol]
    end
    
    methods       
        function Ppc = GasPseudoCritPres(~,dg)
%dg     Gas specific gravity            [ ]
%Ppc    Gas pseudocritical pressure     [kPa]
            Ppc = (46.66 + 1.03*dg - 2.58*dg^2)*1e2;
        end
        
        function Tpc = GasPseudoCritTemp(~,dg)
%dg     Gas specific gravity                    [ ]
%Tpc    Gas pseudocritical temperature          [K]
            Tpc = 93.3 + 180.6*dg - 6.94*dg^2;
        end

        function Zg = EvalGasCompFactor(bo,P,T,dg,Zg)    
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%dg     Gas specific gravity            [ ]
%Ppc    Gas pseudocritical pressure     [kPa]
%Tpc    Gas pseudocritical temperature  [K]
%Zg     Gas compressibility factor      [ ]
            Ppc = GasPseudoCritPres(bo,dg);
            Tpc = GasPseudoCritTemp(bo,dg);
            Zg = 1 + (0.3265 - 1.0700/(T/Tpc) - 0.5339/(T/Tpc)^3 + 0.01569/(T/Tpc)^4 ...
                - 0.05165/(T/Tpc)^5)*(0.27*(P/Ppc)/(Zg*(T/Tpc))) ...
                + (0.5474 - 0.7361/(T/Tpc) + 0.1844/(T/Tpc)^2)*(0.27*(P/Ppc)/(Zg*(T/Tpc)))^2 ...
                - 0.1056*(-0.7361/(T/Tpc) + 0.1844/(T/Tpc)^2)*(0.27*(P/Ppc)/(Zg*(T/Tpc)))^5 ...
                + 0.6134*exp(-0.7210*((0.27*(P/Ppc)/(Zg*(T/Tpc))))^2)* ...
                (1 + 0.7210*(0.27*(P/Ppc)/(Zg*(T/Tpc)))^2)*((0.27*(P/Ppc)/(Zg*(T/Tpc)))^2/(T/Tpc)^3);
        end

        function Zg = GasCompFactor(bo,P,T,dg)
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%dg     Gas specific gravity            [ ]
%Zg     Gas compressibility factor      [ ]
            Zg=1; iter=0; maxiter=100; tol=1e-8;
            while iter<maxiter
                iter=iter+1;
                Zg0=Zg;
                Zg=EvalGasCompFactor(bo,P,T,dg,Zg);
                dev=abs((Zg-Zg0)/Zg0);
                if dev<tol
                    break
                end
            end
        end
        
        function dZgdP = dGasCompFactordP(bo,P,T,dg)
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%dg     Gas specific gravity            [ ]
%Zg     Gas compressibility factor      [ ]
%dZgdP  Zg derivative to P              [1/kPa]
            Ppc = GasPseudoCritPres(bo,dg);
            Tpc = GasPseudoCritTemp(bo,dg);
            Pr=P/Ppc; Tr=T/Tpc;
            Zg = GasCompFactor(bo,P,T,dg);
            dendZg=(1+(1/Zg)*(...
             +(0.3265-1.0700/Tr-0.5339/Tr^3+0.01569/Tr^4-0.05165/Tr^5)*(0.27*Pr/(Zg*Tr))...
             +(0.5474-0.7361/Tr+0.1844/Tr^2)*2*(0.27*Pr/(Zg*Tr))*(0.27*Pr/(Zg*Tr))...
             -0.1056*(-0.7361/Tr+0.1844/Tr^2)*5*(0.27*Pr/(Zg*Tr))^4*(0.27*Pr/(Zg*Tr))...
             +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*(1+0.7210*(0.27*Pr/(Zg*Tr))^2)*((0.27*Pr/(Zg*Tr))^2/Tr^3)*(-0.7210)*2*((0.27*Pr/(Zg*Tr)))*(0.27*Pr/(Zg*Tr))...
             +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*((0.27*Pr/(Zg*Tr))^2/Tr^3)*0.7210*2*((0.27*Pr/(Zg*Tr)))*(0.27*Pr/(Zg*Tr))...
             +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*(1+0.7210*(0.27*Pr/(Zg*Tr))^2)*2/Tr^3*((0.27*Pr/(Zg*Tr)))*(0.27*Pr/(Zg*Tr))));
            dZgdPr=(1/dendZg)*(...
             (0.3265-1.0700/Tr-0.5339/Tr^3+0.01569/Tr^4-0.05165/Tr^5)*(0.27/(Zg*Tr))...
             +(0.5474-0.7361/Tr+0.1844/Tr^2)*2*(0.27*Pr/(Zg*Tr))*(0.27/(Zg*Tr))...
             -0.1056*(-0.7361/Tr+0.1844/Tr^2)*5*(0.27*Pr/(Zg*Tr))^4*(0.27/(Zg*Tr))...
             +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*(1+0.7210*(0.27*Pr/(Zg*Tr))^2)*((0.27*Pr/(Zg*Tr))^2/Tr^3)*(-0.7210)*2*((0.27*Pr/(Zg*Tr)))*(0.27/(Zg*Tr))...
             +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*((0.27*Pr/(Zg*Tr))^2/Tr^3)*0.7210*2*((0.27*Pr/(Zg*Tr)))*(0.27/(Zg*Tr))...
             +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*(1+0.7210*(0.27*Pr/(Zg*Tr))^2)*2/Tr^3*((0.27*Pr/(Zg*Tr)))*(0.27/(Zg*Tr)));
            dZgdP=dZgdPr/Ppc;
        end        

        function dZgdT = dGasCompFactordT(bo,P,T,dg)
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%dg     Gas specific gravity            [ ]
%Zg     Gas compressibility factor      [ ]
%dZgdT  Zg derivative to T              [1/K]
            Ppc = GasPseudoCritPres(bo,dg);
            Tpc = GasPseudoCritTemp(bo,dg);
            Pr=P/Ppc; Tr=T/Tpc;
            Zg = GasCompFactor(bo,P,T,dg);
            dendZg=(1+(1/Zg)*(...
             +(0.3265-1.0700/Tr-0.5339/Tr^3+0.01569/Tr^4-0.05165/Tr^5)*(0.27*Pr/(Zg*Tr))...
             +(0.5474-0.7361/Tr+0.1844/Tr^2)*2*(0.27*Pr/(Zg*Tr))*(0.27*Pr/(Zg*Tr))...
             -0.1056*(-0.7361/Tr+0.1844/Tr^2)*5*(0.27*Pr/(Zg*Tr))^4*(0.27*Pr/(Zg*Tr))...
             +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*(1+0.7210*(0.27*Pr/(Zg*Tr))^2)*((0.27*Pr/(Zg*Tr))^2/Tr^3)*(-0.7210)*2*((0.27*Pr/(Zg*Tr)))*(0.27*Pr/(Zg*Tr))...
             +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*((0.27*Pr/(Zg*Tr))^2/Tr^3)*0.7210*2*((0.27*Pr/(Zg*Tr)))*(0.27*Pr/(Zg*Tr))...
             +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*(1+0.7210*(0.27*Pr/(Zg*Tr))^2)*2/Tr^3*((0.27*Pr/(Zg*Tr)))*(0.27*Pr/(Zg*Tr))));
            dZgdTr=(1/dendZg)*(...
             (-1/Tr)*(...
                (0.3265-1.0700/Tr-0.5339/Tr^3+0.01569/Tr^4-0.05165/Tr^5)*(0.27*Pr/(Zg*Tr))/(-Tr)...
                +(0.5474-0.7361/Tr+0.1844/Tr^2)*2*(0.27*Pr/(Zg*Tr))*(0.27*Pr/(Zg*Tr))/(-Tr)...
                -0.1056*(-0.7361/Tr+0.1844/Tr^2)*5*(0.27*Pr/(Zg*Tr))^4*(0.27*Pr/(Zg*Tr))/(-Tr)...
                +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*(1+0.7210*(0.27*Pr/(Zg*Tr))^2)*((0.27*Pr/(Zg*Tr))^2/Tr^3)*(-0.7210)*2*((0.27*Pr/(Zg*Tr)))*(0.27*Pr/(Zg*Tr))/(-Tr)...
                +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*((0.27*Pr/(Zg*Tr))^2/Tr^3)*0.7210*2*((0.27*Pr/(Zg*Tr)))*(0.27*Pr/(Zg*Tr))/(-Tr)...
                +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*(1+0.7210*(0.27*Pr/(Zg*Tr))^2)*2/Tr^3*((0.27*Pr/(Zg*Tr)))*(0.27*Pr/(Zg*Tr))/(-Tr))...
             +(...
                +(-1.0700/Tr-3*0.5339/Tr^3+4*0.01569/Tr^4-5*0.05165/Tr^5)/(-Tr)*(0.27*Pr/(Zg*Tr))...
                +(-0.7361/Tr+2*0.1844/Tr^2)/(-Tr)*(0.27*Pr/(Zg*Tr))^2 ...
                -0.1056*(-0.7361/Tr+2*0.1844/Tr^2)/(-Tr)*(0.27*Pr/(Zg*Tr))^5 ...
                +0.6134*exp(-0.7210*((0.27*Pr/(Zg*Tr)))^2)*(1+0.7210*(0.27*Pr/(Zg*Tr))^2)*((0.27*Pr/(Zg*Tr))^2/Tr^3)*(-3/Tr))...
             );
            dZgdT=dZgdTr/Tpc;
        end          
        
        function rho_Wstd = WaterStdMassDens(~,ws)
%ws         Water salt mass fraction            [ ]
%rho_Wstd   Water mass density at std. cond.    [kg/m^3]
            rho_Wstd = (62.638 + 0.438603*ws + 1.60074e-3*ws^2)*(0.453592/0.0283168);
        end
        
        function Bo = OilFVF(bo,P,T,d60,dg)
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%d60    Oil specific gravity 60 F       [ ]
%dg     Gas specific gravity            [ ]
%Rs     Solution gas-oil ratio          [Sm^3/Sm^3]
%Bo     Oil formation volume factor     [m^3/Sm^3]
            Rs = SolGasOilRatio(bo,P,T,d60,dg);
           	Bo = 1 + 4.677e-4*(Rs*(0.1589872949/0.0283168)) + 1.751e-5*((T*1.8)-(459.67+60))*(141.5/d60-131.5)/...
                    (dg*(1+5.912e-5*(141.5/d60-131.5)*(40*1.8+32)*log(14.7/114.7))) ...
                - 1.811e-8*(Rs*(0.1589872949/0.0283168))*((T*1.8)-(459.67+60))*(141.5/d60-131.5)/...
                    (dg*(1+5.912e-5*(141.5/d60-131.5)*(40*1.8+32)*log(14.7/114.7))); 
        end

        function dBodP = dOilFVFdP(bo,P,T,d60,dg)
%d60    Oil specific gravity 60 F       [ ]
%dg     Gas specific gravity            [ ]
%Rs     Solution gas-oil ratio          [Sm^3/Sm^3]
%dBodT  Oil FVF derivative to P     	[m^3/Sm^3 kPa]
			dBodRs = dOilFVFdSolGOR(bo,T,d60,dg);
            dRsdP = dSolGasOilRatiodP(bo,P,T,d60,dg);
            dBodP = dBodRs*dRsdP;
        end        
        
        function dBodT = dOilFVFdT(bo,P,T,d60,dg)
%d60    Oil specific gravity 60 F       [ ]
%dg     Gas specific gravity            [ ]
%Rs     Solution gas-oil ratio          [Sm^3/Sm^3]
%dBodT  Oil FVF derivative to T     	[m^3/Sm^3 K]
			Rs=SolGasOilRatio(bo,P,T,d60,dg);
            dBodT = 1.8*(1.751e-5 - 1.811e-8*(Rs*(0.1589872949/0.0283168)))*(141.5/d60-131.5)/ ...
				(dg*(1+5.912e-5*(141.5/d60-131.5)*(40*1.8+32)*log(14.7/114.7)));
        end

        function dBodRs = dOilFVFdSolGOR(~,T,d60,dg)
%T      Temperature                     [K]
%d60    Oil specific gravity 60 F       [ ]
%dg     Gas specific gravity            [ ]
%dBodRs Oil FVF derivative to Sol. GOR  [m^3/Sm^3 K]
			dBodRs = (0.1589872949/0.0283168)*(4.677e-4-1.811e-8*...
                ((T*1.8)-(459.67+60))*(141.5/d60-131.5)/...
                (dg*(1+5.912e-5*(141.5/d60-131.5)*(40*1.8+32)*log(14.7/114.7))));
        end        
        
        function Bg = GasFVF(bo,P,T,dg)
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%dg     Gas specific gravity            [ ]
%Zg     Gas compressibility factor      [ ]
%Bg     Gas formation volume factor     [m^3/Sm^3]
            Zg = GasCompFactor(bo,P,T,dg);
            Bg = (bo.Pref/P)*(T/bo.Tref)*Zg;
        end

        function dBgdP = dGasFVFdP(bo,P,T,dg)
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%dZgdP  Zg derivative to P              [1/kPa]
%Bg     Gas formation volume factor     [m^3/Sm^3]
            Bg = GasFVF(bo,P,T,dg);
            dZgdP = dGasCompFactordP(bo,P,T,dg);
            dBgdP = -Bg/P + (bo.Pref/P)*(T/bo.Tref)*dZgdP;
        end        
        
        function dBgdT = dGasFVFdT(bo,P,T,dg)
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%dZgdT  Zg derivative to P              [1/K]
%Bg     Gas formation volume factor     [m^3/Sm^3]
            Bg = GasFVF(bo,P,T,dg);
            dZgdT = dGasCompFactordT(bo,P,T,dg);      
            dBgdT = Bg/T + (bo.Pref/P)*(T/bo.Tref)*dZgdT;
        end                
        
        function Bw = WaterFVF(~,P,T)
%P      Pressure                            [kPa]
%T      Temperature                         [K]
%DVWP   Water spec. vol. pressure effect    [ ]
%DVWT   Water spec. vol. temperature effect [ ]
%Bw     Water formation volume factor       [m^3/Sm^3]
            DVWP = -(3.58922e-7 + 1.95301e-9*(T*1.8-459.67))*(P*(14.696/101.325)) ...
                -(2.25341e-10 + 1.72834e-13*(T*1.8-459.67))*(P*(14.696/101.325))^2;
            DVWT = -1.0001e-2 + 1.33391e-4*(T*1.8-459.67) + 5.50654e-7*(T*1.8-459.67)^2;
            Bw = (1+DVWP)/(1+DVWT);
        end

        function dBwdP = dWaterFVFdP(~,P,T)
%P      Pressure                            [kPa]
%T      Temperature                         [K]
%DVWP   Water spec. vol. pressure effect    [ ]
%DVWT   Water spec. vol. temperature effect [ ]
%Bw     Water formation volume factor       [m^3/Sm^3]
%dBwdP  Water FVF derivative to P       	[m^3/Sm^3 kPa]
            DVWT = -1.0001e-2 + 1.33391e-4*(T*1.8-459.67) + 5.50654e-7*(T*1.8-459.67)^2;
            dDVWPdP = -(14.696/101.325)*(...
                (3.58922e-7 + 1.95301e-9*(T*1.8-459.67))...
                + 2*(P*14.696/101.325)*(2.25341e-10 + 1.72834e-13*(T*1.8-459.67)));
            dBwdP = dDVWPdP/(1+DVWT);
        end        
        
        function dBwdT = dWaterFVFdT(bo,P,T)
%P      Pressure                            [kPa]
%T      Temperature                         [K]
%DVWP   Water spec. vol. pressure effect    [ ]
%DVWT   Water spec. vol. temperature effect [ ]
%Bw     Water formation volume factor       [m^3/Sm^3]
%dBwdT  Water FVF derivative to T       	[m^3/Sm^3 K]
            DVWT = -1.0001e-2 + 1.33391e-4*(T*1.8-459.67) + 5.50654e-7*(T*1.8-459.67)^2;
			dDVWPdT = -(1.95301e-9)*(P*(14.696/101.325)) - (1.72834e-13)*(P*(14.696/101.325))^2;
            dDVWTdT = 1.33391e-4 + 2*5.50654e-7*(T*1.8-459.67);
			Bw = WaterFVF(bo,P,T);
			dBwdT = 1.8*(-Bw*dDVWTdT + dDVWPdT)/(1+DVWT);
        end
        
        function Rs = SolGasOilRatio(~,P,T,d60,dg)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%d60    Oil specific gravity 60 F   [ ]
%dg     Gas specific gravity        [ ]
%Rs     Solution gas-oil ratio      [Sm^3/Sm^3]
            Rs = ...
                (0.0283168/0.1589872949)*0.0362*(dg*(1 + 5.912e-5*(141.5/d60-131.5)*(40*1.8+32)*log(14.7/114.7)))...
                *(P*(14.696/101.325))^1.0937*exp(25.724*(141.5/d60-131.5)/(T*1.8));
        end

        function dRsdP = dSolGasOilRatiodP(bo,P,T,d60,dg)
%P      Pressure                            [kPa]
%T      Temperature                         [K]
%d60    Oil specific gravity 60 F           [ ]
%dg     Gas specific gravity                [ ]
%Rs     Solution gas-oil ratio              [Sm^3/Sm^3]
%dRsdP  Sol. gas-oil ratio derivative to P  [Sm^3/Sm^3 kPa]
            Rs = SolGasOilRatio(bo,P,T,d60,dg);
            dRsdP = 1.0937*Rs/P;
        end        
        
        function dRsdT = dSolGasOilRatiodT(bo,P,T,d60,dg)
%P      Pressure                            [kPa]
%T      Temperature                         [K]
%d60    Oil specific gravity 60 F           [ ]
%dg     Gas specific gravity                [ ]
%Rs     Solution gas-oil ratio              [Sm^3/Sm^3]
%dRsdT  Sol. gas-oil ratio derivative to T  [Sm^3/Sm^3 K]
            Rs = SolGasOilRatio(bo,P,T,d60,dg);
            dRsdT = -1.8*25.724*(141.5/d60-131.5)/(T*1.8)^2*Rs;
        end
        
        function cPo = OilSpecHeatCap(bo,d60,Tbn)
%d60    Oil specific gravity 60 F   [ ]
%Tbn    Oil normal boiling point    [K]
%cPo    Oil specific heat capacity  [kJ/kg K]
            cPo = (1.05506/0.453592)*1.8*((-1.17126+(0.023722+0.024907*d60)*((Tbn*1.8)^(1/3)/d60) ...
                + (1.14982-0.046535*((Tbn*1.8)^(1/3)/d60))/d60) ...
                + (1/2)*1e-4*(1+0.82463*((Tbn*1.8)^(1/3)/d60))*(1.12172-0.27634/d60)...
                *(707.67^2-(bo.Tref*1.8)^2)/(707.67-(bo.Tref*1.8))...
                - (1/3)*1e-8*(1+0.82463*((Tbn*1.8)^(1/3)/d60))*(2.9027-0.70958/d60)...
                *(707.67^3-(bo.Tref*1.8)^3)/(707.67-(bo.Tref*1.8)));
        end

        function eO = OilSpecEnth(bo,T,d60,Tbn)
%T      Temperature                 [K]
%d60    Oil specific gravity 60 F   [ ]
%Tbn    Oil normal boiling point    [K]
%eO     Oil specific enthalpy       [kJ/kg]
            eO = OilSpecHeatCap(bo,d60,Tbn)*(T-bo.Tref);
        end

        function cPg = GasSpecHeatCap(bo,dg)
%dg     Gas specific gravity        [ ]
%cPg    Gas specific heat capacity  [kJ/kg K]
            cPg = bo.R/(dg*bo.MM_Air)*(...
                (0.7457*dg^2 - 1.3405*dg + 2.222) + ...
                1e-3*(19.894*dg - 1.902)/2*(473.15^2-(bo.Tref)^2)/(473.15-(bo.Tref)) +...
                1e-6*(-6.6562*dg + 1.5103)/3*(473.15^3-(bo.Tref)^3)/(473.15-(bo.Tref))...
                );
        end
        
        function DeOG = GasLatHeat(bo,d60,dg)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%d60    Oil specific gravity 60 F   [ ]
%dg     Gas specific gravity        [ ]
%DeOg   Gas latent heat             [kJ/kg]
            DeOG = bo.R/(dg*bo.MM_Air)*(1/1.8)*(141.5/d60-131.5)*(25.724/1.0937);
        end

        function eGres = GasSpecResEnth(bo,P,T,dg)
%P      Pressure                                [kPa]
%T      Temperature                             [K]
%dg     Gas specific gravity                    [ ]
%Ppc    Gas pseudocritical pressure             [kPa]
%Tpc    Gas pseudocritical temperature          [K]
%eGres  Gas specific residual enthalpy          [kJ/kg]
            Ppc = GasPseudoCritPres(bo,dg);
            Tpc = GasPseudoCritTemp(bo,dg);
            eGres = bo.R/(dg*bo.MM_Air)*Tpc*(...
                (P/Ppc)*((0.083-1.097/(T/Tpc)^1.6) + (0.1432*dg - 0.0669)*(0.139-0.894/(T/Tpc)^4.2)) - ...
                (bo.Pref/Ppc)*((0.083-1.097/(bo.Tref/Tpc)^1.6) + (0.1432*dg - 0.0669)*(0.139-0.894/(bo.Tref/Tpc)^4.2)) ...
                );
        end

        function deGresdP = dGasSpecResEnthdP(bo,T,dg)
%P      	Pressure                                [kPa]
%T      	Temperature                             [K]
%dg     	Gas specific gravity                    [ ]
%Ppc    	Gas pseudocritical pressure             [kPa]
%Tpc       	Gas pseudocritical temperature          [K]
%deGresdT  	eGres derivative to T                   [kJ/kg K]
            Ppc = GasPseudoCritPres(bo,dg);
            Tpc = GasPseudoCritTemp(bo,dg);
            deGresdP = bo.R/(dg*bo.MM_Air)*Tpc*(...
                (1/Ppc)*((0.083-1.097/(T/Tpc)^1.6) + (0.1432*dg - 0.0669)*(0.139-0.894/(T/Tpc)^4.2)));
        end

        function deGresdT = dGasSpecResEnthdT(bo,P,T,dg)
%P      	Pressure                                [kPa]
%T      	Temperature                             [K]
%dg     	Gas specific gravity                    [ ]
%Ppc    	Gas pseudocritical pressure             [kPa]
%Tpc       	Gas pseudocritical temperature          [K]
%deGresdT  	eGres derivative to T                   [kJ/kg K]
            Ppc = GasPseudoCritPres(bo,dg);
            Tpc = GasPseudoCritTemp(bo,dg);
			deGresdT = bo.R/(dg*bo.MM_Air)*(...
                (P/Ppc)*(-(1.6)*(-1.097/(T/Tpc)^2.6)-(4.2)*(0.1432*dg-0.0669)*(-0.894/(T/Tpc)^5.2)));	
        end
            
        function eGig = IdealGasSpecEnth(bo,T,dg)
%T      Temperature                     [K]
%dg     Gas specific gravity            [ ]
%eGig   Ideal gas specific enthalpy     [kJ/kg]
            eGig = GasSpecHeatCap(bo,dg)*(T-bo.Tref);
        end

        function cPw = WaterSpecHeatCap(bo,P,T,ws)
%P          Pressure                                                [kPa]
%T          Temperature                                             [T]
%ws         Water salt mass fraction                                [ ]
%rho_Wstd   Water mass density at std. cond.                        [kg/m^3]
%DVWP       Water spec. vol. pressure effect                        [ ]
%IDVWT      Integral of water spec. vol. temperature effect         [ ]
%IDVWTSQ    Integral of water spec. vol. temperature effect square  [ ]
%cPw        Water specific heat capacity                            [kJ/kg K]
            rho_Wstd = WaterStdMassDens(bo,ws);
            DVWP = -(3.58922e-7 + 1.95301e-9*(T*1.8-459.67))*(P*(14.696/101.325)) ...
                -(2.25341e-10 + 1.72834e-13*(T*1.8-459.67))*(P*(14.696/101.325))^2;
            IDVWT = -1.0001e-2 + 1.33391e-4/2*(707.67^2-(bo.Tref*1.8)^2)/(707.67-(bo.Tref*1.8)) ...
                + 5.50654e-7/3*(707.67^3-(bo.Tref*1.8)^3)/(707.67-(bo.Tref*1.8));
            IDVWTSQ = 1.0001^2e-4 - 1.33391e-6*(707.67^2-(bo.Tref*1.8)^2)/(707.67-(bo.Tref*1.8)) ...
                + 2.25966e-9*(707.67^3-(bo.Tref*1.8)^3)/(707.67-(bo.Tref*1.8)) ...
                + 3.67261e-11*(707.67^4-(bo.Tref*1.8)^4)/(707.67-(bo.Tref*1.8)) ...
                + 6.06440e-14*(707.67^5-(bo.Tref*1.8)^5)/(707.67-(bo.Tref*1.8));
            cPw = (4.18 - ws*(rho_Wstd)/(1+DVWP)*(4.396e-3 - 4.8e-6*ws*(rho_Wstd)/(1+DVWP))) ...
                - ws*(rho_Wstd)/(1+DVWP)*(4.396e-3 - 2*4.8e-6*ws*(rho_Wstd)/(1+DVWP))*IDVWT ...
                + 4.8e-6*(ws*(rho_Wstd)/(1+DVWP))^2*IDVWTSQ;
        end

        function eW = WaterSpecEnth(bo,P,T,ws)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%ws     Water salt mass fraction    [ ]
%eW     Water specific enthalpy     [kJ/kg]
            eW = WaterSpecHeatCap(bo,P,T,ws)*(T-bo.Tref);
        end
        
    end
end