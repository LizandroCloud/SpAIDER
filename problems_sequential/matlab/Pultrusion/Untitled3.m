function [T, C] =Pultrude(Line_Speed, Thickness, Charge_Temperature,T1, T2, T3, T4, Die_Length);
tcondy=10.4;
sp_ht=1010;
sp_htg=1970
densityg=2450;
densityr=1350
glass_VF=0.6;
h=500;
density=glass_vf*(densityg)+(1+glass_vf)*densityr;
mass_f=glass_vf*desityg/density;
cp=mass_f*sp_htg+(1-mass_f)*sp_ht;
A1=1.044E13;
E1_k=124;
A2=3.686E8;
E2_k=73.435;
m=0.67;
n=1.37;
Q=210;
glass_vf=0.6;
min_Fo=0.1;
min_DT=2;
min_DC=0.1;
E1=E1_k*1000;
E2=E2_k*1000;
rQ=Q*(1-0.60*2450/densityg);
step_size=5;
NumSteps=round(Die_Length/step_size);
OneThird=round(NumSteps/3);
FirstThird=OneThird;
TwoThirds=2*OneThird;
First_Slope=(T2-T1)/OneThird;
Second_Slope=(T3-T2)/OneThird;
Third_Slope=(T4-T3)/OneThird;
Speed=Line_Speed*1000/60;
alpha_p=tcondy/density/sp_ht;
time=step_size/Speed;
delta_thickness=sqrt(alpha_p*time/min_Fo*1E6);
inner=Thickness/delta_thickness;
slices=round(inner);
delta_thickness=Thickness/(slices);
fo=alpha_p*time/(delta_thickness*delta_thickness)*1e6;
Biot=h*delta_thickness/1000/tcondy;
rQ=Q*(1-glass_VF*2450/density);
for i=1:NumSteps
for j=1:slices
T(i,j)=Charge_Temperature;
C(i,j)=0.0;
end
end
for i=1:NumSteps
if i<=OneThird Die_Temperature=T1+First_Slope*i;
end
if ((i>OneThird) & (i<=TwoThirds))
Die_Temperature=T2+Second_Slope*(i-OneThird);
end
if i>TwoThirds Die_Temperature = T3+ Third_Slope*(i-TwoThirds);
end
T(i+1,1)=2*Biot*fo*Die_Temperature+ 2*fo*T(i,2)+(1-2*fo*(Biot+1))*T(i,1);
for j=2:slices-1
    T(i+1,j)= T(i,j-1)*fo + T(i,j+1)*fo+ T(i,j)*(1-2*fo);
end
T(i+1,slices)=2*Biot*fo*Die_Temperature+2*fo*T(i,slices-1)+(1-2*fo*(Biot+1))*T(i,slices);
j=0;
while j<slices
j=j+1;
This_Temperature=T(i,j);
This_Conversion=C(i,j);
New_Conversion=RK4(This_Conversion,This_Temperature,time);
Reaction_T=(New_Conversion-This_Conversion)*rQ*1000/sp_ht;%
if New_Conversion< This_Conversion
New_Conversion=This_Conversion;
end
%
T(i+1,j)=T(i+1,j)+Reaction_T;
C(i+1,j)=New_Conversion;
end
end
plot (T(:,1))
hold on
xlabel('Die Lengthx5 mm')
ylabel('Temperature C')
title('Temp vs Die Length')
axis([0 180 0 180])
hold on
plot(C(:,1))
xlabel('Die Lengthx5 mm')
ylabel('Conversion')
title('Conversion vs Die Length')
end