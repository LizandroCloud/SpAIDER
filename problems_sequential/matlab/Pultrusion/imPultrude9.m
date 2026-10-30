function[T, C, Die_Temperature]= imPultrude9(Line_Speed, Thickness,Charge_Temperature, Die_Length);
close all;
fftest=input('enter value :');
tcondy=10.4;
sp_ht=1010;
sp_htg=1970
densityg=2450;
densityr=1350
glass_vf=0.6;
h=500;
density=glass_vf*(densityg)+(1+glass_vf)*densityr;
mass_f=glass_vf*densityg/density;
cp=mass_f*sp_htg+(1-mass_f)*sp_ht;
A1=1.044E13;
E1_k=124;
A2=3.686E8;
E2_k=73.435;
m=0.67;
n=1.37;
Area=0.018;
dx=0.01;
tconds=17;
sdensity=7850;
sp_hts=450;
Q=1150;
glass_vf=0.6;
min_Fo=0.1;
min_DT=2;
min_DC=0.1;
E1=E1_k*1000;
E2=E2_k*1000;
rQ=Q*(1-0.65*2450/density);
step_size=5;
NumSteps=round(Die_Length/step_size);
Speed=Line_Speed*1000/60;
alpha_p=tcondy/(density*sp_ht);
alpha_s=tconds/(sdensity*sp_hts);
time=step_size/Speed;
delta_thickness=sqrt(alpha_p*time/min_Fo*1E6);
inner=Thickness/delta_thickness;
slices=round(inner);
delta_thickness=Thickness/(slices);
fo=alpha_p*time/(delta_thickness*delta_thickness)*1e6;
Biot=h*delta_thickness/1000/tcondy;
rQ=Q*(1-glass_vf*2450/density);
for i=1:fftest
for j=1:8
T(i,j)=Charge_Temperature;
C(i,j)=0.0;
end
end
Die_Temperature(1)=180;
Heat(1)=1150;
p(1)=1;
A=[(1+2*fo+2*fo*Biot),-2*fo,0,0,0,0,0,0;-fo,(1+2*fo),- ...
fo,0,0,0,0,0;0,-fo,(1+2*fo),-fo,0,0,0,0;0,0,-fo,(1+2*fo),- ...
fo,0,0,0;0,0,0,-fo,(1+2*fo),-fo,0,0;0,0,0,0,-fo,(1+2*fo),- ...
fo,0;0,0,0,0,0,-fo,(1+2*fo),-fo;0,0,0,0,0,0,- ...
2*fo,(1+2*fo+2*fo*Biot)];
for i=1:fftest-1
if i==1
B=[(T(i,1)+(2*fo*Biot*Die_Temperature(i)));T(i,2);T(i,3);T(i,4);T(i,5);
    T(i,6);T(i,7);(T(i,8)+(2*fo*Biot*Die_Temperature(i)))];
D=inv(A)*B;
Dans= D';
end
if i>1
for j=1:8
for b=1:10
    if i== b*NumSteps
B=[Charge_Temperature;Charge_Temperature;Charge_Temperature;Charge_T
emperature;Charge_Temperature;Charge_Temperature;Charge_Temperature;
Charge_Temperature];
D=B;
else
B=D;
end
end
end
B(1,1)=B(1,1)+(2*fo*Biot*Die_Temperature(i));
B(8,1)=B(8,1)+(2*fo*Biot*Die_Temperature(i));
D=inv(A)*B;
Dans=[Dans ; D'];
end
Die_Temperature(i+1)=(Heat(i)/(sdensity*sp_hts*Area*dx))- ...
(alpha_s*(abs(Die_Temperature(i)-(T(i))))/dx^2)+Die_Temperature(i);
if Die_Temperature(i+1)>180
Heat(i+1)=0;
p(i+1)=0;
else
Heat(i+1)=1150;
p(i+1)=1;
end
j=0;
while j<8
j=j+1;
This_Temperature=Dans(i,j);
This_Conversion=C(i,j);
New_Conversion=RK4(This_Conversion,This_Temperature,time);
Reaction_T=(New_Conversion- ...
This_Conversion)*rQ*1000/sp_ht;%
if New_Conversion< This_Conversion
New_Conversion=This_Conversion;
end
T(i+1,j)=Dans(i,j)+Reaction_T;
C(i+1,j)=New_Conversion;
end
plot(Die_Temperature(1,:))
end
% bar(p(1,:));
 plot (T)
% hold on
% xlabel('Die Length mm')
ylabel('Temperature C')
% title('Temp vs Die Length')
% hold on
 plot(C(:,1))
 xlabel('Die Length mm')
 ylabel('Conversion')
end