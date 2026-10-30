function Y=RK4(Conversion, Temperature, Step);
S1=Step*KamalRate(Conversion, Temperature);
S2=Step*KamalRate(Conversion+S1/2, Temperature);
S3=Step*KamalRate(Conversion+S2/2, Temperature);
S4=Step*KamalRate(Conversion+S3, Temperature);
Y=Conversion+(S1+2*S2+2*S3+S4)/6;