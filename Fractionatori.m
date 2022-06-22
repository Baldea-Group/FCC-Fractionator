function [Temperatureout,Vaporout,xcout,Liqout,Holdout,EnthalLout,EnthalVout,LN,HN,LCO,ELC2,ETC4,ETC5,ETC6,Ttrack1,Ttrack2,yout,Valvesf]=Fractionatori(xfra,ufra,xc,MV,SP,products,errord,yin,Xfil,dist,cont)
%**************************************************************************************************************************************************************************************************************************************************************************************************************************
%Fractionator dynamic model (from top to bottom)
%Inputs:
%xfra-->Vector of states (The vector is concatenated as a column: [Hold-up; Enthalpy liquid; Enthalpy vapor; Liquid flowrate].  Each element of the subvectors (e.g  Hold-up, Enthalpy liquid) represents a fractionator stage. The elements of the subvectors are ordered from the top stage (1) to the bottom stage (20))
%ufra-->Vector of decision variables for NLP (The vector is concatenated as a column: [Vapor flowrate; Temperature]. Each element of the subvectors (e.g  Temperature, Vapor flowrate) represents a fractionator stage. The elements of the subvectors are ordered from the top stage (1) to the bottom stage (20))
%xc-->Matrix of liquid composition (Each row (i) of the matrix represents the liquid composition of the stage. Each column (j) represents the component composition. The rows are ordered from the top stage (1) to the bottom stage (20). The columns are ordered from lightest component (Methane) to heaviest one (VGO)) 
%MV-->Vector of manipulated variables [reflux,flowrate of water at condenser,duties for PA, flowrates for PA]
%SP-->Vector of set points for controllers [accum. level (Lfra),overhead temp. (Tfra), cut point HN (Thnt) cut point LCO (Tlcot)]
%products-->Vector of fractionator liquid products flowrate [Distillate,HN,LCO]
%errord-->Vector of controllers integral error [accum. level, overhead temp., cut point HN, cut point LCO]
%yin-->Matrix of vapor composition (Each row (i) of the matrix represents the liquid composition of the stage. Each column (j) represents the component composition. The rows are ordered from the top stage (1) to the bottom stage (20). The columns are ordered from lightest component (Methane) to heaviest one (VGO)) 
%Xfil-->Vector of filtered FCC reactor outputs [Feed flowrate to fractionator, Temperature of feed stream, Composition of feed stream,Pressure of feed stream]
%dist-->Vector of disturbances
%cont-->Time step
%Outputs:
%Temperatureout-->Vector of temperature of stages. The elements of the vector are ordered from the top stage (1) to the bottom stage (20)
%Vaporout-->Vector of vapor flowrate of stages. The elements of the vector are ordered from the top stage (1) to the bottom stage (20)
%xcout-->Matrix of liquid composition (Each row (i) of the matrix represents the liquid composition of the stage. Each column (j) represents the component composition. The rows are ordered from the top stage (1) to the bottom stage (20). The columns are ordered from lightest component (Methane) to heaviest one (VGO)) 
%Liqout-->Vector of liquid flowrate of stages. The elements of the vector are ordered from the top stage (1) to the bottom stage (20)
%Holdout-->Vector of liquid hold-up of stages. The elements of the vector are ordered from the top stage (1) to the bottom stage (20)
%EnthalLout-->Vector of liquid Enthalpy of stages. The elements of the vector are ordered from the top stage (1) to the bottom stage (20)
%EnthalVout-->Vector of vapor Enthalpy of stages. The elements of the vector are ordered from the top stage (1) to the bottom stage (20)
%LN-->LN flowrate (molar)
%HN-->HN flowrate (molar)
%LCO-->LCO flowrate (molar)
%ELC2-->Error of controller LC2
%ETC4-->Error of controller TC4
%ETC5-->Error of controller TC5
%ETC6-->Error of controller TC6
%Ttrack1-->98% Heavy tail cut point for HN
%Ttrack2-->98% Heavy tail cut point for LCO
%yout-->Matrix of vapor composition (Each row (i) of the matrix represents the liquid composition of the stage. Each column (j) represents the component composition. The rows are ordered from the top stage (1) to the bottom stage (20). The columns are ordered from lightest component (Methane) to heaviest one (VGO)) 
%Valvesf-->Vector of valve positions [V9,V8,V10,V11]
%**************************************************************************************************************************************************************************************************************************************************************************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%Tower/Solver Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=20; %Number of stages
k=10;%Number of components
delta=10/3600;% Integration step (hours)
options = optimoptions('fsolve','FunctionTolerance',1e-6,'OptimalityTolerance',1e-6,'Display','none');
Mwin=[85.1347950005991;120.941794467086;206.951432493898;292.185529267655;378.894679684485;446.094358949534];

%%%%%%%%%%%%%%%%%%%%%Definition of feed conditions%%%%%%%%%%%%%%%%%%%%%%%%
Feedr=Xfil(1);%Kmol/h 
Tfeed=Xfil(2);%K 
xfeed=[Xfil(3);Xfil(4);Xfil(5);Xfil(6);Xfil(7);Xfil(8);Xfil(9);Xfil(10);Xfil(11);Xfil(12)]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%Definition of MV's%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reflux=MV(1);
effi=dist(4);%Condenser efficiency
CpCoolant=4.18;%(KJ/KgK)
Tcoolant=(dist(1)-32)*(5/9)+273.15;%(K)
UDutycond=(effi*(MV(2)*CpCoolant)*(ufra(N+2)-Tcoolant))/1000; 
UDuty=[UDutycond;MV(3);MV(4)];  %Duties from top to bottom    (KMJ/h)
UPump=[MV(5);MV(6)]; %2 pump arounds top to bottom

%%%%%%%%%%%%%%%%%%%%%%%%%%%Transfrom SP to K%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SP(2)=(SP(2)-32)*(5/9)+273.15;
SP(3)=(SP(3)-32)*(5/9)+273.15;
SP(4)=(SP(4)-32)*(5/9)+273.15;

%%%%%%%%%%%%%%%%%%%%%%Definition of products%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Distillate=products(1);
Fhn=products(2);
Flco=products(3);

%%%%%%%%%%%%%%%%%%%%%Definition of NLP variables%%%%%%%%%%%%%%%%%%%%%%%%%%%
UVap=ufra(1:N);
UTemp=ufra(N+1:N+N);
Ux=xc; %Matrix [20,10] format

%%%%%%%%%%%%%%%%%%%%%%%%Definition of states%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UMl=xfra(1:N);
EnL=xfra(N+1:N+N);
EnV=xfra(N+N+1:N+N+N);
Liq=xfra(N+N+N+1:N+N+N+N);

%%%%%%%%%%%%%%%%%%Pressure computation per stage%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pcond=Xfil(13)-(5.48/14.503773773);%Pressure condenser
for i=1:N-2
  Pa(i)=Xfil(13)+(0.1*0.0689476*(i));   
end
Pfra=[Pcond,Xfil(13),Pa];
Pfeed=Xfil(13)+(5.4/14.503773773);

%%%%%%%%%%%%%%%%%%%%%%%STAGE MASS AND ENERGY BALANCES%%%%%%%%%%%%%%%%%%%%%%
Temperatureout=[];
Vaporout=[];
Liqout=[];
Holdout=[];
EnthalLout=[];
EnthalVout=[];
xcout=[];
yout=[];
for l=1:N

if l==1    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Condenser%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
j=l; 
Feed=UVap(l+1);
z=[];
for i=1:k
a=(UVap(l+1)*yin(l+1,i))/Feed;
z=[z;a];
end

H=(UVap(l+1)*EnV(l+1))/Feed;        
[Us, fval]= fsolve(@stage,[UTemp(l);UVap(l);Ux(l,:)'],options,Feed,z,H,delta,UMl(l),Pfra(l),Ux(l,:),EnL(l),j,Distillate,UDuty,UPump,Liq(1));
[Mhold,liquid,y,HL,HV]=holds(Us,Feed,delta,UMl(l),Pfra(l),j,Distillate,UPump,Liq(1));    

%%%%%%%%%%%%%%%%%%%%%Temperature Controller (eTC4)%%%%%%%%%%%%%%%%%%%%%%%%%
eTC4=-1*(SP(2)-UTemp(l+1));
 
KcTC4=4;                            
TaoTC4=0.6; 
V8nom=50;
Rliquidbase=810.630566926004;

if cont>=1  
V8a=V8nom+KcTC4*(eTC4+ errord(2)/TaoTC4);
V8b=max(5,V8a); 
V8=min(95,V8b); 
else
eTC4=0;    
V8=V8nom;
end 

ETC4=eTC4*delta+errord(2);
liquid=Rliquidbase*(V8/V8nom);

%%%%%%%%%%%%%%%%%%%%%%%%%%Level Controller (eLC2)%%%%%%%%%%%%%%%%%%%%%%%%%%
eLC2=(SP(1)-UMl(l))*-1;
 
KcLC2=8;                                
TaoLC2=1;
V9nom=50;
Distillatebase=810.630566926004/Reflux;

if cont>=1  
 V9a=V9nom+KcLC2*(eLC2+errord(1)/TaoLC2);
V9b=max(5,V9a);
V9=min(95,V9b); 
else
 eLC2=0;    
 V9=V9nom;    
end 

LN=Distillatebase*((V9/V9nom))^(0.5); 
ELC2=eLC2*delta+errord(1);

elseif l==6    
%%%%%%%%%%%%%%%%%%%%%%%%%%Side HN & Controller%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
j=l; 
Feed=UVap(l+1)+Liqout(end);
z=[];
for i=1:k
a=(UVap(l+1)*yin(l+1,i)+Liqout(end)*xcout(end,i))/Feed;
z=[z;a];
end

H=(UVap(l+1)*EnV(l+1)+Liqout(end)*EnthalLout(end))/Feed;
[Us, fval]= fsolve(@stage,[UTemp(l);UVap(l);Ux(l,:)'],options,Feed,z,H,delta,UMl(l),Pfra(l),Ux(l,:),EnL(l),j,Fhn,UDuty,UPump,Flco);
[Mhold,liquid,y,HL,HV]=holds(Us,Feed,delta,UMl(l),Pfra(l),j,Fhn,UPump,Flco); 

%%%%%%%%%%%%%%%%%%%%%%%%%%HN controller (eTC5)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xshc=max(0,(0.02-Ux(l,8)));
Nxg=Ux(l,8)+xshc;
Ttrack1=(6.135762930946595*1.0e+02*Ux(l,8)+5.303668836827089*1.0e+02*xshc)/Nxg;

%alpha=(100/200)
eTC5=(SP(3)-Ttrack1)*1;
 
KcTC5=5;                        
TaoTC5=1;
p10nom=100;
HNbase=101.874173498128;

if cont>=1 
 p10a=p10nom+KcTC5*(eTC5+errord(3)/TaoTC5);
 p10b=max(10,p10a);
 p10=min(190,p10b);  
else
 eTC5=0;    
 p10=p10nom;    
end 

ETC5=eTC5*delta+errord(3);
HN=HNbase*(p10/p10nom)^(0.5); 
V10=p10*(1/2);

elseif l==9 
%%%%%%%%%%%%%%%%%%%%%%%%%%%PA Top Column%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=l; 
Feed=UVap(l+1)+Liqout(end)+UPump(1);
z=[];
for i=1:k
a=(UVap(l+1)*yin(l+1,i)+Liqout(end)*xcout(end,i)+UPump(1)*yin(l+2,i))/Feed;
z=[z;a];
end

H=(UVap(l+1)*EnV(l+1)+Liqout(end)*EnthalLout(end)+(UPump(1)*EnV(l+2)-UDuty(2)))/Feed;
[Us, fval]= fsolve(@stage,[UTemp(l);UVap(l);Ux(l,:)'],options,Feed,z,H,delta,UMl(l),Pfra(l),Ux(l,:),EnL(l),j,Distillate,UDuty,UPump,Flco);
[Mhold,liquid,y,HL,HV]=holds(Us,Feed,delta,UMl(l),Pfra(l),j,Distillate,UPump,Flco); 

elseif l==13  
%%%%%%%%%%%%%%%%%%%%Side LCO & Controller%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
j=l; 
Feed=UVap(l+1)+Liqout(end);
z=[];
for i=1:k
a=(UVap(l+1)*yin(l+1,i)+Liqout(end)*xcout(end,i))/Feed;
z=[z;a];
end

H=(UVap(l+1)*EnV(l+1)+Liqout(end)*EnthalLout(end))/Feed;
[Us, fval]= fsolve(@stage,[UTemp(l);UVap(l);Ux(l,:)'],options,Feed,z,H,delta,UMl(l),Pfra(l),Ux(l,:),EnL(l),j,Distillate,UDuty,UPump,Flco);
[Mhold,liquid,y,HL,HV]=holds(Us,Feed,delta,UMl(l),Pfra(l),j,Distillate,UPump,Flco); 

%%%%%%%%%%%%%%%%%%%%%%%%%%LCO controller (eTC6)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsdc=max(0,(0.02-Ux(l,9)));
Nxd=Ux(l,9)+xsdc;
Ttrack2=(6.833520046817131*1.0e+02*Ux(l,9)+6.135762930946595*1.0e+02*xsdc)/Nxd;

%alpha=(100/200)
eTC6=(SP(4)-Ttrack2)*1;  
 
KcTC6=2;                                      
TaoTC6=1;           
p11nom=100;
Lcobase=168.410080092446;

    if cont>=1 
     p11a=p11nom+KcTC6*(eTC6+errord(4)/TaoTC6);
     p11b=max(10,p11a);
     p11=min(190,p11b);  
    else
     eTC6=0;    
     p11=p11nom;    
    end 

ETC6=eTC6*delta+errord(4);
LCO=Lcobase*(p11/p11nom)^(0.5);
V11=p11*(1/2);

elseif l==18
%%%%%%%%%%%%%%%%%%%%%%PA Bottom Column%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=l; 
Feed=UVap(l+1)+Liqout(end)+UPump(end);
z=[];
for i=1:k
a=(UVap(l+1)*yin(l+1,i)+Liqout(end)*xcout(end,i)+UPump(end)*yin(l+2,i))/Feed;
z=[z;a];
end

H=(UVap(l+1)*EnV(l+1)+Liqout(end)*EnthalLout(end)+(UPump(end)*EnV(l+2)-UDuty(end)))/Feed;
[Us, fval]= fsolve(@stage,[UTemp(l);UVap(l);Ux(l,:)'],options,Feed,z,H,delta,UMl(l),Pfra(l),Ux(l,:),EnL(l),j,Distillate,UDuty,UPump,Flco);
[Mhold,liquid,y,HL,HV]=holds(Us,Feed,delta,UMl(l),Pfra(l),j,Distillate,UPump,Flco);     

elseif l==20
%%%%%%%%%%%%%%%%%%%%%%%%%%Bottom Stage%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[HLini,HVini]=EnthalpyB(Tfeed,Pfeed,xfeed',Mwin);   
j=l;    

Feed=Feedr+Liqout(end);
z=[];
for i=1:k
a=(Feedr*xfeed(i)+Liqout(end)*xcout(end,i))/Feed;
z=[z;a];
end

H=(Feedr*HVini+Liqout(end)*EnthalLout(end))/Feed;
[Us, fval]= fsolve(@stage,[UTemp(l);UVap(l);Ux(l,:)'],options,Feed,z,H,delta,UMl(l),Pfra(l),Ux(l,:),EnL(l),j,Distillate,UDuty,UPump,Flco);
[Mhold,liquid,y,HL,HV]=holds(Us,Feed,delta,UMl(l),Pfra(l),j,Distillate,UPump,Flco);

else
j=l;    
Feed=UVap(l+1)+Liqout(end);
z=[];
for i=1:k
a=(UVap(l+1)*yin(l+1,i)+Liqout(end)*xcout(end,i))/Feed;
z=[z;a];
end

H=(UVap(l+1)*EnV(l+1)+Liqout(end)*EnthalLout(end))/Feed;
[Us, fval]= fsolve(@stage,[UTemp(l);UVap(l);Ux(l,:)'],options,Feed,z,H,delta,UMl(l),Pfra(l),Ux(l,:),EnL(l),j,Distillate,UDuty,UPump,Flco);
[Mhold,liquid,y,HL,HV]=holds(Us,Feed,delta,UMl(l),Pfra(l),j,Distillate,UPump,Flco);   
 
end

Temperatureout=[Temperatureout;Us(1)];
Vaporout=[Vaporout;Us(2)];
Liqout=[Liqout;liquid];
Holdout=[Holdout;Mhold];
EnthalLout=[EnthalLout;HL];
EnthalVout=[EnthalVout;HV];
xcout=[xcout;Us(3:end)'];
yout=[yout;y'];

end

Valvesf=[V9,V8,V10,V11];
end

function er=stage(U,Feed,z,H,delta,Mhold0,P,x,H0,j,UF1,UDuty,UPump,UF2)
%NLP. Solves the mass and energy balance equations 
%Inputs:
%U-->Optimization variables [Temperature, Vapor and liquid composition of stage at future time step]
%Feed->Feed to the stage (considered an adiabatic flash)
%z-->Composition of feed stream (considered an adiabatic flash)
%H-->Enthalpy of feed stream (considered an adiabatic flash)
%delta-->Integration step
%Mhold0-->Hold-up of the stage at current time step
%P-->Pressure of the stage
%x-->Liquid composition of stage at current time step
%H0-->Liquid enthalpy of stage at current time step
%j-->Stage #
%UF1-->Flow of product (either LN(distillate),HN or LCO)
%UDuty-->Duty of PA
%UPump-->Flow of PA
%UF2-->Flow of product (either LN(distillate),HN or LCO)
%Outputs:
%er-->Residuals to be minimized

%%%%%%%%%%%%%%%%%%%%%Definition of Optimization variables%%%%%%%%%%%%%%%%%%
UTemp=U(1); %Temperature 
UVap=U(2); %Vapor
Ux=U(3:end)'; %Composition
k=10;%Number of components
Mwin=[85.1347950005991;120.941794467086;206.951432493898;292.185529267655;378.894679684485;446.094358949534];

%%%%%%%%%%%%%%%%%%%%%%%%Thermodynamic computations%%%%%%%%%%%%%%%%%%%%%%%%%
[y,Zl,Keq,HL,HV]=Enthalpy(UTemp,P,Ux,Mwin);

%%%%%%%%%%%%%%%%%%%%%%%Definition of HTC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The HTC for the accumulator is used to compute its hold-up value from
%the SS solution. For the dynamic operation, the hold-up value is 
%determined by the level and temperature controllers 
%(eLC2 & eTC4 respectively)
if j==1
tao=5*60/3600;    
elseif j==20
tao=1.2*60/3600;    
else
tao=1*60/3600;    
end

if j==1 %%%%%%%%%%%%%%%%%%%%%%Accumulator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The hold-up is determined by the controllers (eLC2 & eTC4) 
%Hold up Balance
Mhold=(Feed-UF1-UF2-UVap+(Mhold0/delta))*(delta);

Mass=[];
%Mass balance
for i=1:k
mas=(Feed*(z(i)-Ux(i))-UVap*(y(i)-Ux(i)))/(Mhold)-((Ux(i)-x(i))/(delta));
Mass=[Mass; mas];
end

%Energy Balance
Energy=(Feed*(H-HL)-UVap*(HV-HL)-UDuty(1))/(Mhold)-((HL-H0)/(delta));

%Mole Balance
moles=((ones(1,k)*y)-(Ux*ones(k,1)));

%Error Vector
er=[Mass;Energy;moles]; 

elseif j==6 %%%%%%%%%%%%%%%%%%%%%%HN side%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Hold up computation
Mhold=(Feed-UF1-UVap+(Mhold0/delta))*(tao*delta/(delta+tao));

Mass=[];
%Mass balance
for i=1:k
mas=(Feed*(z(i)-Ux(i))-UVap*(y(i)-Ux(i)))/(Mhold)-((Ux(i)-x(i))/(delta));
Mass=[Mass; mas];
end

%Energy Balance
Energy=(Feed*(H-HL)-UVap*(HV-HL))/(Mhold)-((HL-H0)/(delta));

%Mole Balance
moles=((ones(1,k)*y)-(Ux*ones(k,1)));
 
%Error Vector
er=[Mass;Energy;moles]; 

elseif j==11 %%%%%%%%%%%%%%%%%%%%%%%Top PA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hold up computation
Mhold=(Feed-UPump(1)-UVap+(Mhold0/delta))*(tao*delta/(delta+tao));

Mass=[];
%Mass balance
for i=1:k
mas=(Feed*(z(i)-Ux(i))-UVap*(y(i)-Ux(i))-UPump(1)*(y(i)-Ux(i)))/(Mhold)-((Ux(i)-x(i))/(delta));
Mass=[Mass; mas];
end

%Energy Balance
Energy=(Feed*(H-HL)-UVap*(HV-HL)-UPump(1)*(HV-HL))/(Mhold)-((HL-H0)/(delta));

%Mole Balance
moles=((ones(1,k)*y)-(Ux*ones(k,1)));

%Error Vector
er=[Mass;Energy;moles]; 

elseif j==13 %%%%%%%%%%%%%%%%%%%%%%%%LCO side%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hold up computation
Mhold=(Feed-UF2-UVap+(Mhold0/delta))*(tao*delta/(delta+tao));

Mass=[];
%Mass balance
for i=1:k
mas=(Feed*(z(i)-Ux(i))-UVap*(y(i)-Ux(i)))/(Mhold)-((Ux(i)-x(i))/(delta));
Mass=[Mass; mas];
end

%Energy Balance
Energy=(Feed*(H-HL)-UVap*(HV-HL))/(Mhold)-((HL-H0)/(delta));

%Mole Balance
moles=((ones(1,k)*y)-(Ux*ones(k,1)));

%Error Vector
er=[Mass;Energy;moles]; 

elseif j==20 %%%%%%%%%%%%%%%%%%%%%%%%%Bottom PA%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hold up computation
Mhold=(Feed-UPump(end)-UVap+(Mhold0/delta))*(tao*delta/(delta+tao));

Mass=[];
%Mass balance
for i=1:k
mas=(Feed*(z(i)-Ux(i))-UVap*(y(i)-Ux(i))-UPump(end)*(y(i)-Ux(i)))/(Mhold)-((Ux(i)-x(i))/(delta));
Mass=[Mass; mas];
end

%Energy Balance
Energy=(Feed*(H-HL)-UVap*(HV-HL)-UPump(end)*(HV-HL))/(Mhold)-((HL-H0)/(delta));

%Mole Balance
moles=((ones(1,k)*y)-(Ux*ones(k,1)));

%Error Vector
er=[Mass;Energy;moles]; 

else %%%%%%%%%%%%%%%%%%%%%%%%%Other stages%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hold up computation
Mhold=(Feed-UVap+(Mhold0/delta))*(tao*delta/(delta+tao));

Mass=[];
%Mass balance
for i=1:k
mas=(Feed*(z(i)-Ux(i))-UVap*(y(i)-Ux(i)))/(Mhold)-((Ux(i)-x(i))/(delta));
Mass=[Mass; mas];
end

%Energy Balance
Energy=(Feed*(H-HL)-UVap*(HV-HL))/(Mhold)-((HL-H0)/(delta));

%Moles
 moles=((ones(1,k)*y)-(Ux*ones(k,1)));
 er=[Mass;Energy;moles]; 

end
end

function [Mhold,liquid,y,HL,HV]=holds(U,Feed,delta,Mhold0,P,j,UF1,UPump,UF2)
%Computes the hold-up, liquid, vapor composition, and the enthalpies of the two phases 
%Inputs:
%U-->Solution of NLP [Temperature, Vapor and liquid composition of stage at future time step]
%Feed->Feed to the stage (considered an adiabatic flash)
%delta-->Integration step
%Mhold0-->Hold-up of the stage at current time step
%P-->Pressure of the stage
%j-->Stage #
%UF1-->Flow of product (either LN(distillate),HN or LCO)
%UPump-->Flow of PA
%UF2-->Flow of product (either LN(distillate),HN or LCO)
%Outputs:
%Mhold-->Hold-up of the stage at future time step
%liquid-->Liquid flowrate of the stage at future time step
%y-->Vapor composition of the stage at future time step
%HL-->Liquid Enthalpy of the stage at future time step
%HV-->Vapor Enthalpy of the stage at future time step

%%%%%%%%%%%%%%%%%%%%%%%%Definition of Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UTemp=U(1); %Temperature 
UVap=U(2); %Vapor
Ux=U(3:end)'; %Composition
Mwin=[85.1347950005991;120.941794467086;206.951432493898;292.185529267655;378.894679684485;446.094358949534];

%%%%%%%%%%%%%%%%%%%%%%%%Thermodynamic computations%%%%%%%%%%%%%%%%%%%%%%%%%
[y,Zl,Keq,HL,HV]=Enthalpy(UTemp,P,Ux,Mwin);

%%%%%%%%%%%%%%%%%%%%%%%Definition of HTC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The HTC for the accumulator is used to compute its hold-up value from
%the SS solution. For the dynamic operation, the hold-up value is 
%determined by the level and temperature controllers 
%(eLC2 & eTC4 respectively)
if j==1
tao=5*60/3600;    
elseif j==20
tao=1.2*60/3600;    
else
tao=1*60/3600;    
end

if j==1 %%%%%%%%%%%%%%%%%%%%%%Accumulator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The hold-up is determined by the controllers (eLC2 & eTC4)
%Hold up computation
Mhold=(Feed-UF1-UF2-UVap+(Mhold0/delta))*(delta);
liquid=UF2;

elseif j==6 %%%%%%%%%%%%%%%%%%%%%%HN side%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hold up computation
Mhold=(Feed-UF1-UVap+(Mhold0/delta))*(tao*delta/(delta+tao));
liquid=(Mhold)/(tao);

elseif j==11 %%%%%%%%%%%%%%%%%%%%%%%Top PA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hold up computation
Mhold=(Feed-UPump(1)-UVap+(Mhold0/delta))*(tao*delta/(delta+tao));
liquid=(Mhold)/(tao);

elseif j==13 %%%%%%%%%%%%%%%%%%%%%%%%LCO side%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hold up computation
Mhold=(Feed-UF2-UVap+(Mhold0/delta))*(tao*delta/(delta+tao));
liquid=(Mhold)/(tao);  

elseif j==20 %%%%%%%%%%%%%%%%%%%%%%%%%Bottom PA%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hold up computation
Mhold=(Feed-UPump(end)-UVap+(Mhold0/delta))*(tao*delta/(delta+tao)); 
liquid=(Mhold)/(tao);

else %%%%%%%%%%%%%%%%%%%%%%%%%Other stages%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hold up computation
Mhold=(Feed-UVap+(Mhold0/delta))*(tao*delta/(delta+tao));
liquid=(Mhold)/(tao);
end

end