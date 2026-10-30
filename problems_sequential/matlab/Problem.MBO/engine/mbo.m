classdef mbo
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

        function Zg = EvalGasCompFactor(mbo,P,T,dg,Zg)    
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%dg     Gas specific gravity            [ ]
%Ppc    Gas pseudocritical pressure     [kPa]
%Tpc    Gas pseudocritical temperature  [K]
%Zg     Gas compressibility factor      [ ]
            Ppc = GasPseudoCritPres(mbo,dg);
            Tpc = GasPseudoCritTemp(mbo,dg);
            Zg = 1 + (0.3265 - 1.0700/(T/Tpc) - 0.5339/(T/Tpc)^3 + 0.01569/(T/Tpc)^4 ...
                - 0.05165/(T/Tpc)^5)*(0.27*(P/Ppc)/(Zg*(T/Tpc))) ...
                + (0.5474 - 0.7361/(T/Tpc) + 0.1844/(T/Tpc)^2)*(0.27*(P/Ppc)/(Zg*(T/Tpc)))^2 ...
                - 0.1056*(-0.7361/(T/Tpc) + 0.1844/(T/Tpc)^2)*(0.27*(P/Ppc)/(Zg*(T/Tpc)))^5 ...
                + 0.6134*exp(-0.7210*((0.27*(P/Ppc)/(Zg*(T/Tpc))))^2)* ...
                (1 + 0.7210*(0.27*(P/Ppc)/(Zg*(T/Tpc)))^2)*((0.27*(P/Ppc)/(Zg*(T/Tpc)))^2/(T/Tpc)^3);
        end
        
        function rho_Wstd = WaterStdMassDens(~,ws)
%ws         Water salt mass fraction            [ ]
%rho_Wstd   Water mass density at std. cond.    [kg/m^3]
            rho_Wstd = (62.638 + 0.438603*ws + 1.60074e-3*ws^2)*(0.453592/0.0283168);
        end
        
        function Bo = OilFVF(mbo,P,T,d60,dg)
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%d60    Oil specific gravity 60 F       [ ]
%dg     Gas specific gravity            [ ]
%Rs     Solution gas-oil ratio          [Sm^3/Sm^3]
%dgs    Solution gas specific gravity   [ ]
%Bo     Oil formation volume factor     [m^3/Sm^3]
            Rs = SolGasOilRatio(mbo,P,T,d60,dg);
            dgs = SolGasSpecGrav(mbo,P,T,dg);
           	Bo = 1 + 4.677e-4*(Rs*(0.1589872949/0.0283168)) + 1.751e-5*((T*1.8)-(459.67+60))*(141.5/d60-131.5)/...
                    (dgs*(1+5.912e-5*(141.5/d60-131.5)*(40*1.8+32)*log(14.7/114.7))) ...
                - 1.811e-8*(Rs*(0.1589872949/0.0283168))*((T*1.8)-(459.67+60))*(141.5/d60-131.5)/...
                    (dgs*(1+5.912e-5*(141.5/d60-131.5)*(40*1.8+32)*log(14.7/114.7))); 
        end

        function dBodT = dOilFVFdT(~,Rs,d60,dgs)
%d60    Oil specific gravity 60 F       [ ]
%dgs    Solution gas specific gravity   [ ]
%Rs     Solution gas-oil ratio          [Sm^3/Sm^3]
%dBodT  Oil FVF derivative to T     	[m^3/Sm^3 K]
			dBodT = 1.8*(1.751e-5 - 1.811e-8*(Rs*(0.1589872949/0.0283168)))*(141.5/d60-131.5)/ ...
				(dgs*(1+5.912e-5*(141.5/d60-131.5)*(40*1.8+32)*log(14.7/114.7)));
        end

%        function Bg = GasFVF(mbo,P,T,dg)
%P      Pressure                        [kPa]
%T      Temperature                     [K]
%dg     Gas specific gravity            [ ]
%Zg     Gas compressibility factor      [ ]
%Bg     Gas formation volume factor     [m^3/Sm^3]
%            Zg = GasCompFactor(mbo,P,T,dg);
%            Bg = (mbo.Pref/P)*(T/mbo.Tref)*Zg;
%        end

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

        function dBwdT = dWaterFVFdT(mbo,P,T)
%P      Pressure                            [kPa]
%T      Temperature                         [K]
%DVWP   Water spec. vol. pressure effect    [ ]
%DVWT   Water spec. vol. temperature effect [ ]
%Bw     Water formation volume factor       [m^3/Sm^3]
%dBwdT  Water FVF derivative to T       	[m^3/Sm^3 K]
            DVWT = -1.0001e-2 + 1.33391e-4*(T*1.8-459.67) + 5.50654e-7*(T*1.8-459.67)^2;
			dDVWPdT = -(1.95301e-9)*(P*(14.696/101.325)) - (1.72834e-13)*(P*(14.696/101.325))^2;
            dDVWTdT = 1.33391e-4 + 2*5.50654e-7*(T*1.8-459.67);
			Bw = WaterFVF(mbo,P,T);
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

        function Rv = VapOilGasRatio(~,P,T,d60,dg,Ycg)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%d60    Oil specific gravity 60 F   [ ]
%dg     Gas specific gravity        [ ]
%Ycg    Condensate-Gas yield        [Sm^3/SMm^3]
%Rv     Vaporized oil-gas ratio     [Sm^3/SMm^3]
            Rv = Ycg*(4.7e-7*(P*(14.696/101.325))^2 - 1.4e-3*(P*(14.696/101.325)) + 2.3)* ...
                exp(-5.438e-3*dg)*exp(-6.8e8*d60/(T*1.8)^3);
        end

        function d60G = HypOilSpecGrav(mbo,P,T,d60)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%d60    Oil specific gravity 60 F   [ ]
%d60G   Hyp. oil spec. grav. 60 F   [ ]
            d60G = d60*(mbo.Pref/P)^0.010*exp(100/450 - 100/T);
        end

        function TbnG = HypOilNormBolPnt(mbo,P,T,d60,Tbn)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%d60    Oil specific gravity 60 F   [ ]
%Tbn    Oil normal boiling point    [K]
%d60G   Hyp. oil spec. grav. 60 F   [ ]
            d60G = HypOilSpecGrav(mbo,P,T,d60);
            TbnG = Tbn*(d60G/d60)^3;
        end
        
        function dgs = SolGasSpecGrav(mbo,P,T,dg)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%dg     Gas specific gravity        [ ]
%dgs    Solution gas spec. grav.    [ ]
            dgs = dg/((mbo.Pref/P)^0.005*exp(600/450 - 600/T));
        end

        function Ycgs = SolCondGasRatio(mbo,P,T,Ycg)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%Ycg    Condensate-Gas yield        [Sm^3/SMm^3]
%Ycgs   Sol. condensate-gas yield   [Sm^3/SMm^3]
            Ycgs = Ycg/((mbo.Pref/P)^0.005*exp(600/450 - 600/T));
        end
        
        function cPo = OilSpecHeatCap(mbo,d60,Tbn)
%d60    Oil specific gravity 60 F   [ ]
%Tbn    Oil normal boiling point    [K]
%cPo    Oil specific heat capacity  [kJ/kg K]
            cPo = (1.05506/0.453592)*1.8*((-1.17126+(0.023722+0.024907*d60)*((Tbn*1.8)^(1/3)/d60) ...
                + (1.14982-0.046535*((Tbn*1.8)^(1/3)/d60))/d60) ...
                + (1/2)*1e-4*(1+0.82463*((Tbn*1.8)^(1/3)/d60))*(1.12172-0.27634/d60)...
                *(707.67^2-(mbo.Tref*1.8)^2)/(707.67-(mbo.Tref*1.8))...
                - (1/3)*1e-8*(1+0.82463*((Tbn*1.8)^(1/3)/d60))*(2.9027-0.70958/d60)...
                *(707.67^3-(mbo.Tref*1.8)^3)/(707.67-(mbo.Tref*1.8)));
        end

        function eO = OilSpecEnth(mbo,T,d60,Tbn)
%T      Temperature                 [K]
%d60    Oil specific gravity 60 F   [ ]
%Tbn    Oil normal boiling point    [K]
%eO     Oil specific enthalpy       [kJ/kg]
            eO = OilSpecHeatCap(mbo,d60,Tbn)*(T-mbo.Tref);
        end

        function cPg = GasSpecHeatCap(mbo,dg)
%dg     Gas specific gravity        [ ]
%cPg    Gas specific heat capacity  [kJ/kg K]
            cPg = mbo.R/(dg*29)*(...
                (0.7457*dg^2 - 1.3405*dg + 2.222) + ...
                1e-3*(19.894*dg - 1.902)/2*(473.15^2-(mbo.Tref)^2)/(473.15-(mbo.Tref)) +...
                1e-6*(-6.6562*dg + 1.5103)/3*(473.15^3-(mbo.Tref)^3)/(473.15-(mbo.Tref))...
                );
        end
        
        function DeOG = GasLatHeat(mbo,P,T,d60,dg)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%d60    Oil specific gravity 60 F   [ ]
%dg     Gas specific gravity        [ ]
%d60G   Hyp. oil spec. grav. 60 F   [ ]
%DeOg   Gas latent heat             [kJ/kg]
            d60G = HypOilSpecGrav(mbo,P,T,d60);
            DeOG = mbo.R/(dg*29)*(1/1.8)*(141.5/d60G-131.5)*(25.724/1.0937);
        end

        function eGres = GasSpecResEnth(mbo,P,T,dg)
%P      Pressure                                [kPa]
%T      Temperature                             [K]
%dg     Gas specific gravity                    [ ]
%Ppc    Gas pseudocritical pressure             [kPa]
%Tpc    Gas pseudocritical temperature          [K]
%eGres  Gas specific residual enthalpy          [kJ/kg]
            Ppc = GasPseudoCritPres(mbo,dg);
            Tpc = GasPseudoCritTemp(mbo,dg);
            eGres = mbo.R/(dg*29)*Tpc*(...
                (P/Ppc)*((0.083-1.097/(T/Tpc)^1.6) + (0.1432*dg - 0.0669)*(0.139-0.894/(T/Tpc)^4.2)) - ...
                (mbo.Pref/Ppc)*((0.083-1.097/(mbo.Tref/Tpc)^1.6) + (0.1432*dg - 0.0669)*(0.139-0.894/(mbo.Tref/Tpc)^4.2)) ...
                );
        end
            
        function eGig = IdealGasSpecEnth(mbo,T,dg)
%T      Temperature                     [K]
%dg     Gas specific gravity            [ ]
%eGig   Ideal gas specific enthalpy     [kJ/kg]
            eGig = GasSpecHeatCap(mbo,dg)*(T-mbo.Tref);
        end

        function cPw = WaterSpecHeatCap(mbo,P,T,ws)
%P          Pressure                                                [kPa]
%T          Temperature                                             [T]
%ws         Water salt mass fraction                                [ ]
%rho_Wstd   Water mass density at std. cond.                        [kg/m^3]
%DVWP       Water spec. vol. pressure effect                        [ ]
%IDVWT      Integral of water spec. vol. temperature effect         [ ]
%IDVWTSQ    Integral of water spec. vol. temperature effect square  [ ]
%cPw        Water specific heat capacity                            [kJ/kg K]
            rho_Wstd = WaterStdMassDens(mbo,ws);
            DVWP = -(3.58922e-7 + 1.95301e-9*(T*1.8-459.67))*(P*(14.696/101.325)) ...
                -(2.25341e-10 + 1.72834e-13*(T*1.8-459.67))*(P*(14.696/101.325))^2;
            IDVWT = -1.0001e-2 + 1.33391e-4/2*(707.67^2-(mbo.Tref*1.8)^2)/(707.67-(mbo.Tref*1.8)) ...
                + 5.50654e-7/3*(707.67^3-(mbo.Tref*1.8)^3)/(707.67-(mbo.Tref*1.8));
            IDVWTSQ = 1.0001^2e-4 - 1.33391e-6*(707.67^2-(mbo.Tref*1.8)^2)/(707.67-(mbo.Tref*1.8)) ...
                + 2.25966e-9*(707.67^3-(mbo.Tref*1.8)^3)/(707.67-(mbo.Tref*1.8)) ...
                + 3.67261e-11*(707.67^4-(mbo.Tref*1.8)^4)/(707.67-(mbo.Tref*1.8)) ...
                + 6.06440e-14*(707.67^5-(mbo.Tref*1.8)^5)/(707.67-(mbo.Tref*1.8));
            cPw = (4.18 - ws*(rho_Wstd)/(1+DVWP)*(4.396e-3 - 4.8e-6*ws*(rho_Wstd)/(1+DVWP))) ...
                - ws*(rho_Wstd)/(1+DVWP)*(4.396e-3 - 2*4.8e-6*ws*(rho_Wstd)/(1+DVWP))*IDVWT ...
                + 4.8e-6*(ws*(rho_Wstd)/(1+DVWP))^2*IDVWTSQ;
        end

        function eW = WaterSpecEnth(mbo,P,T,ws)
%P      Pressure                    [kPa]
%T      Temperature                 [K]
%ws     Water salt mass fraction    [ ]
%eW     Water specific enthalpy     [kJ/kg]
            eW = WaterSpecHeatCap(mbo,P,T,ws)*(T-mbo.Tref);
        end
        
    end
end