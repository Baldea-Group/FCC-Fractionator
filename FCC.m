function [delta,yp]=FCC(xfcc,dist,ufcc,Flpg,Tcondenser,cont)
%************************************************************************************************************************************************************************
%Model based on McFarlane et al theoretical developments and its Matlab implementation (references [40] and [43] respectively)
%The names and symbols of the original variables/parameters (McFarlane study) IN GENERAL are preserved (some modifications are performed, primarily in the valves naming)
%New variables/parameters are introduced to account for the improved FCC model.
%Inputs:
%xfcc-->Vector of states
%dist-->Vector of disturbances
%ufcc-->Vector of manipulated variables
%Flpg-->Flow of LPG
%Tcondenser-->Condenser temperature
%cont-->Time step
%Outputs:
%delta-->Derivatives
%yp-->Process "measurements"
%************************************************************************************************************************************************************************
global Fair FH vs rhocdl Crgc zbed P6 epse Treg Fsg 
global Tcyc o2cyc cocyc co2cyc
global xo2 xco xco2
global Wris Tr Frgc P4

 %>>>>>>>>>Define distubances

   Tatm = dist(1); %Ambient temperature disturbance
   API = dist(2); %Feed quality disturbance
   T1   = dist(3); %Feed temperature disturbance
   APInominal=25; %Nominal API (APInominal)
   
%'>>>>>>>Feed System Constants' 
   taufb=200.0; 
   DHfu=1000.0*(API/APInominal)^0.5; 
   UAf=25.0;
   a1=0.15;     
   a2=200.0;    
   taufo=60.0;
   F5nom=34; %Nominal Flow of fuel to furnace

 %>>>>>>>>Wet Gas Compressor Constants.'

	k13=0.01;	
	k11=1200.0212842935; 
	Pvru=75.0;	
    
 %>>>>>>>>Lift Air BLower Constants.'

	vslip=2.2;	
	taufil=40.0;	
	hlift=34.0;
	sb=5950.0;	
	samin=5000.0;	
	kavg=1.39;
	etapla=1.0;	
	Tdla=225.0;		
	k9=10.0;
	Vdla=200.0;	
	klift=5.0;		
	k8=5.0;
    
    %>>>>>>>>Combustion Air Blower Constants.'

	Tdca=190.0;
	kca=40.0;		
	k6=250.0;	
	k7=15.0;
	Vcms=200.0;	
	Vdca=1000.0;

 %>>>>>>>>Regenerator Constants.'

	rgarea=590.0; 
	sparea=7.0;		
	cpCO2=11.0;
	cpair=7.08;		
	cpc=0.31;		
	cpCO=7.28;
	cpH2O=8.62;		
	cpN2=7.22;		
	cpO2=7.62; 
	k14=1.3464; 			
	MI=200000.0;	
	Qe=556.0;
	zcyc=45.0;		
	zsp=13.0;		
	R=10.73;					
	zlp=11.0;
	rhoprt=68.0;		
	rhoc=45.0;		
	Patm=14.7;		
	delH1=46368.0;
	fof=424;			
	hsp=20.0;		
	Tair=270.0;
	delH2=169080.0;	
	delhH=60960.0;	
	Tf=459.6;
    
    %>>>>>>Catalyst Recirculation:'

	areaur=3.7;  
	areaus=5.2;     
	Astrp=60.0;   
	Alp=8.73;    
	Lur=56.0;    
	Lus=56.0;       
	Etap=155.0;
	Eloi=124.5;  
	Estrp=130.0;    
	Elift=134.0;   
	furgc=17.0;     
	fusc=47.0;

% '>>>>>>>>Reactor Riser:'

	Mcpeff=10000.0; 
	Tbase=1100.0;   
	Tref=999.0;
	Tbasef=700.0;   
	dTstrp=35.0;    
	Hcoke=0.075;
	cpsv=0.80;      
	cpfv=0.81*(API/APInominal)^(-0.1);      
	cpfl=0.82*(API/APInominal)^(-0.1);   
	Qfr=309.0*(API/APInominal)^(-0.15);      
	Qsr=412.0;      
	rsarea=9.6;   
	hris=60.0;              
	rhov=0.57; 
    
% '>>>>>>>>Reactor Pressure:'

	k12=0.5; 
	dPfrac=9.5;   

%>>>>>>>>Define Variable for Common block
      
	T3=xfcc(1); %Furnace firebox temperature
	T2=xfcc(2); %Temperature of fresh feed entering the reactor (Tpre IN MANUSCRIPT)
	P7=xfcc(3); %WGC suction pressure
	P5=xfcc(4); %Fractionator pressure (Pfra IN MANUSCRIPT)
	P3=xfcc(5); %Lift air blower discharge pressure
	P6=xfcc(6); %Regenerator pressure (Preg IN MANUSCRIPT)
	rho=xfcc(7); %Density of catalyst in lift pipe
	P2=xfcc(8); %CAB discharge pressure
	Csc=xfcc(9); %Weight fraction of coke on spent catalyst
    Crgc=xfcc(10); %Weight fraction of coke on regenerated catalyst
    Treg=xfcc(11); %Temperature of regenerator (Treg IN MANUSCRIPT)
	Wsp=xfcc(12); %Inventory catalyst in regenerator standpipe
	Wreg=xfcc(13); %Inventory catalyst in regenerator
	Rn=xfcc(14); %Quantity of gas
	Wr=xfcc(15); %Inventory catalyst in reactor (Lrea IN MANUSCRIPT)
	Tr=xfcc(16); %Temperature of reactor riser (Trea IN MANUSCRIPT)
    Fair=xfcc(17); %Flow of air into regenerator 
	Pblp=xfcc(18); %Pressure at bottom of lift pipe
	P1=xfcc(19); %CAB suction pressure
    PreHeatE=xfcc(20); %Preheater temperature controller integral error
    MainFracE=xfcc(21); %Fractionator temperature controller integral error
    RegenE=xfcc(22); %Regenerator pressure controller integral error
    RegenTE=xfcc(23); %Regenerator temperature controller integral error
    CatIE=xfcc(24); %Reactor catalyst inventory controller integral error
    ReacTE=xfcc(25); %Reactor temperature controller integral error
    GFpvgo=xfcc(26);  %CST MB for VGO
    GFp1=xfcc(27); %CST MB for PC1
    GFp2=xfcc(28); %CST MB for PC2
    GFp3=xfcc(29); %CST MB for PC3
    GFp4=xfcc(30); %CST MB for PC4
    GFc5=xfcc(31); %CST MB for C5+
    GFb=xfcc(32); %CST MB for Butane
    GFp=xfcc(33); %CST MB for Propane
    GFe=xfcc(34); %CST MB for Ethane
    GFm=xfcc(35); %CST MB for Methane
    Fwg=xfcc(36); %Filter for LPG
    Fcoke=xfcc(37); %Filter for Coke

   GTotal=[GFpvgo;GFp1;GFp2;GFp3;GFp4;GFc5;GFb;GFp;GFe;GFm]; 
   P4		= P5 + dPfrac;
   deltP	= P6 - P4;
   	
 %>>>>>>>>>Define manipulated variables 
	F3=ufcc(1); %Flow of fresh feed
	F4=ufcc(2); %Flow of slurry
	V12=ufcc(3); %Combustion air blower suction valve position
	V13=ufcc(4); %Combustion air blower vent valve position
	V14=ufcc(5); %Lift air blower vent valve position
	V15=ufcc(6); %Spill air valve position
	Vlift=ufcc(7); %Lift air blower steam
	V5=ufcc(8); %Wet gas flare valve position
	V16=ufcc(9); %WGC vent valve position
    SPPH=ufcc(10); %Set point preheater (Temperature)
    SPMF=ufcc(11); %Set point fractionator (Pressure) 
    SPRG=ufcc(12); %Set point regenerator (Pressure) 
    SPRGT=ufcc(13); %Set point regenerator (Temperature) 
    SPCI=ufcc(14); %Set point reactor (Catalyst inventory) 
	SPRT=ufcc(15); %Set point reactor (Riser temperature) 
    
    V13	=max(0.0d0,min(1.0d0,V13));
	V14	=max(0.0d0,min(1.0d0,V14));
	Vlift	=max(0.0d0,min(1.2d0,Vlift));
    
%>>>>>>>>>Nonlinear Valve Macro.'
 
	if(V12<=0.5)
   	V12=0.3*V12;
	else
   	V12=exp(-3.79424*(1.0d0-V12));
	end
  
	if(V13<=0.5)
 	 	V13=0.3*V13;
	else
  		V13=exp(-3.79424*(1.0d0-V13));
	end

%*************************************************************'
%                  TIME DERIVATIVES                            '
%**************************************************************'

%>>>>>>>>Independent variables

Tsc    = Tr - dTstrp;
Prgb   = P6 + Wreg/(144.0*rgarea);

%%%%%%%%%%%%%%%%% Controller Fractionator Pressure (PC2) %%%%%%%%%%%%%%%%%%
%alpha=1, hence, p(t)=valve(t)
ePC2=(SPMF-P5)*-1;    

KcPC2=1500;                                
TaoPC2=200;                                                                         
V4nom=50;

if cont>=1   
V4a=V4nom+KcPC2*(ePC2+MainFracE/TaoPC2);
V4b=max(5,V4a);
V4=min(95,V4b);
else
ePC2=0;    
V4=V4nom; 
end

FV11   = (k11/2)*(V4/V4nom)^2*sqrt(max(0.0d0,P5 - P7)); 
FV12   = k12*V5*sqrt(max(0.0d0,P5 - Patm));
FV13   = k13*V16*Pvru;
dMainFracE=ePC2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Controller Regenerator Temperature (TC3) %%%%%%%%%%%%%%%%
%alpha=(100/160)
eTC3=(SPRGT-Treg)*1;
    
KcTC3=1.5;                    
TaoTC3=800;            
p6nom=50;

if cont>=1
p6a=p6nom+KcTC3*(eTC3+RegenTE/TaoTC3);
p6b=max(8,p6a);
p6=min(152,p6b);
else
eTC3=0;    
p6=p6nom;    
end                

F7     = kca*(p6/p6nom)^1*sqrt(max(0.0d0,P2-Prgb)); %(Fa IN MANUSCRIPT)
V6=p6*(100/160);
dRegenTE=eTC3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FV7    = k7*V13*sqrt(max(0.0d0,P2-Patm));
F10    = k9*V15*sqrt(max(0.0d0,P3-Prgb));
Fwgi    =Flpg; 

%%%%%%%%%%%%%%%%%%%%% Controller Regenerator Pressure (PC1)%%%%%%%%%%%%%%%%
%alpha=1, hence, p(t)=valve(t)
ePC1=(SPRG-P6)*-1;     

KcPC1=1500;                             
TaoPC1=200;                              
V7nom=50;

if cont>=1
V7a=V7nom+KcPC1*(ePC1+RegenE/TaoPC1);
V7b=max(10,V7a);
V7=min(95,V7b);
else
ePC1=0;    
V7=V7nom;    
end 

Fsg    = (1/2)*(V7/V7nom)*k14*sqrt(max(0.0d0,P6-Patm)); %(Ffg IN MANUSCRIPT)
dRegenE=ePC1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhog   = 0.0933*P6/(Treg + Tf);   
F9     = klift*sqrt(max(0.0d0,P3 - Pblp));
vs     = ((Fsg + Fair)/2.0)*(1.0/(rhog*rgarea));
epsf   = 0.332 + 0.06*vs;
rhocdn = rhoprt*(1.0 - epsf);
rhocdl = -0.878 + 0.582*vs;
zbed   = min(zcyc,(2.85 + 0.8*vs + (Wreg - rhocdl*rgarea*...
               zcyc)/(rgarea*rhocdn))*(1.0/(1.0 - rhocdl/rhocdn)));
epse   = max(epsf,epsf + (1.904 + 0.363*vs - 0.048*vs*vs)/zbed);
epse   = min(1.0d0,epse);

%>>>>>>>Spent Catalyst Ubend.'

dPsc   = 144.0*(P4-Pblp) + ((Wr/Astrp) + (Estrp-Elift)*rhoc);

%%%%%%%%%%%%%%%%%% Controller Catalyst Inventory (LC1)%%%%%%%%%%%%%%%%%%%%%                   
%alpha=1, hence, p(t)=valve(t)
eLC1=(SPCI-Wr)*-1;

KcLC1=0.0001;                       
TaoLC1=1000;  
V3nom=50;

if cont>=1
V3a = V3nom+KcLC1*(eLC1+CatIE/TaoLC1);
V3b=max(5,V3a);
V3=min(95,V3b);
else
eLC1=0;    
V3=V3nom;
end

vsc    = (dPsc*areaus)/(Lus*fusc);
Fsc    = (V3/V3nom)^2*vsc*areaus*rhoc; %(Fsc IN MANUSCRIPT)
dCatIE=eLC1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FH     = Fsc*(Csc - Crgc)*Hcoke; 

%>>>>>>>>Preheat System.'

DTin     = T3 - T1;
DTout    = T3 - T2;
Tlm      = (DTin - DTout)/(log(DTin/DTout));

%%%%%%%%%%%%%%%%%%%%%%% Controller Preheater (TC1)%%%%%%%%%%%%%%%%%%%%%%%%%
%alpha=(118.75/100)
eTC1=(SPPH-T2)*1;

KcTC1=5;                          
TaoTC1=100;                        
p1nom=50;

if cont>=1   
p1a=p1nom+KcTC1*(eTC1+PreHeatE/TaoTC1);
p1b=max(20,p1a); 
p1=min(80,p1b); 
else
eTC1=0;    
p1=p1nom;
end    
F5=(p1/p1nom)*F5nom; %(Ff IN MANUSCRIPT)
V1=p1*(118.75/100);
Qloss    = a1*F5*T3 - a2;
T2ss     = T1 + (UAf*Tlm/F3);
dPreHeatE=eTC1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%>>>>>>>>>Wet Gas Compressor.'
	  
Crw      = Pvru/P7;
Hwg      = 6884936.03341729*(Crw^2 - 1.0d0); 
Fsucwg   = 106859.958742798 + sqrt(max(0.0d0,592504703729978-0.15*Hwg^2.0)); 
F11      = 2.636d-6*Fsucwg*P7;

%>>>>>>>>>Lift Air Blower.'

M      = (kavg-1.0)/kavg/etapla;
sa     = samin + 1100.0*Vlift;
xp5p   = max(0.0d0,P3);     
Pbdla  = ((xp5p^M-Patm^M)*(sb/sa)^2.0 + Patm^M)^(1.0d0/M);
Fbsla  = 8600.0 + sqrt(max(0.0d0,2.582d8-1.068d5*Pbdla^2.0));
Fsucla = Fbsla*sa/sb; 
F8     = 0.0451*Patm*Fsucla/(Tatm+Tf);
FV8    = k8*V14*sqrt(max(0.0d0,P3-Patm));
rhoag  = 29.0*P6/R/(Tsc+Tf);
vairl  = F9/Alp/rhoag;
vcatl  = max(vairl-vslip,Fsc/Alp/rhoprt);

%>>>>>>>>>Combustion Air Blower.'

Pbca   = 14.7*P2/P1;
Fsucca = 45000.0+sqrt(max(0.0d0,1.5813d9-1.2491d6*Pbca^2.0));
F6     = 0.045*Fsucca*P1/(Tatm + Tf);
Fv6    = k6*V12*sqrt(Patm-P1);

%>>>>>>>Regenerated Catalyst Ubend.'
%%%%%%%%%%%%%%%%%%% Controller Reactor Temperature (TC2) %%%%%%%%%%%%%%%%%%
%alpha=1, hence, p(t)=valve(t)
eTC2=(SPRT-Tr)*1;

KcTC2=0.5;                       
TaoTC2=200;  
V2nom=50;

if cont>=1
V2a= V2nom+KcTC2*(eTC2+ReacTE/TaoTC2);
V2b=max(5,V2a);
V2=min(95,V2b);    
else
eTC2=0;    
V2=V2nom;
end  

eta1   = (Wsp/sparea) + (Etap-Eloi)*rhoc;
eta2   = areaur/(Lur*furgc);
eta3   = 144.0*eta2*(P6-P4) + eta2*eta1;
rho1   = rhov*rhoprt;
alpha1 = areaur*rhoc*eta3;
alpha2 = areaur*rhoc*eta2*hris;
beta1  = rhoprt*(F3+F4) + alpha2*rho1 - alpha1*rhov;
beta2  = (F3+F4)*(alpha2*rho1 - alpha1*rhoprt);
Frgc   = (V2/V2nom)^2*(-beta1 + ((beta1^2) - 4*rhov*beta2)^(0.5))/2.0/rhov; %(Frgc IN MANUSCRIPT)
dReacTE=eTC2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vris   = (F3 + F4)/rhov + Frgc/rhoprt;
rhoris = (F3 + F4 + Frgc)/vris;
Wris   = (Frgc*rsarea*hris)/vris;
Prb    = P4 + (rhoris*hris)/144.0;

a0p=[2.24167820089635;0;0;0;0;0;0;0;0;0;0];%Initial Concentration of feed to reactor (mol/Kg gas)
[ypfr,taoF]=PFR(a0p,F3,F4,vris,API);%Riser model (PFR subroutine)

 dGFpvgo=taoF*(ypfr(1)-GFpvgo);% CST model
 dGFp1=taoF*(ypfr(2)-GFp1);% CST model
 dGFp2=taoF*(ypfr(3)-GFp2);% CST model
 dGFp3=taoF*(ypfr(4)-GFp3);% CST model
 dGFp4=taoF*(ypfr(5)-GFp4);% CST model
 dGFc5=taoF*(ypfr(6)-GFc5);% CST model
 dGFb=taoF*(ypfr(7)-GFb);% CST model
 dGFp=taoF*(ypfr(8)-GFp);% CST model
 dGFe=taoF*(ypfr(9)-GFe);% CST model
 dGFm=taoF*(ypfr(10)-GFm);% CST model

MWp5=446.094358949534;
MWp4=378.894679684485;
MWp3=292.185529267655;
MWp2=206.951432493898;
MWp1=120.941794467086;
MWc5=85.1347950005991;
MWb=58.12;
MWp=44.1;
MWe=30.07;
MWm=16.04;
MWc=400;

MWT=[MWp5;MWp4;MWp3;MWp2;MWp1;MWc5;MWb;MWp;MWe;MWm;MWc]; %(g/mol)
 
yfinal=([GTotal;ypfr(end)].*MWT)*(F3+F4)*(1/2.20462)*(1/453.592); %Mass of each component (lb/s)

Fcokei=yfinal(end);

%>>>>>>>>>Regenerator Calculations:'

zposition %Regenerator temperature profile and MB

Qair   = Fair*cpair*(Tair - Tbase);
Qh     = FH*delhH;
Qc     = Fair*(cocyc*delH1 + co2cyc*delH2);
Qsc    = Fsc*cpc*(Tsc-Tbase);
Qfg    = (Fair*(o2cyc*cpO2 + cocyc*cpCO + co2cyc*cpCO2 + 0.79*cpN2) + 0.5*FH*cpH2O)*(Tcyc - Tbase);
Qrgc   = Frgc*cpc*(Treg - Tbase);
Qin    = Qair + Qh + Qc + Qsc;
Qout   = Qfg + Qrgc + Qe;
hspp   = hsp-(Wsp/rhoc/sparea);
Fspadj = 20.0*(3.0 - min(3.0d0,hspp));
Fsp    = fof*(sparea^0.5)*(zbed - zsp) - 4925.0 - Fspadj;
splev  = Wsp/(rhoc*sparea);
Vregg  = rgarea*zcyc - rgarea*zbed*(1.0-epse);

dn     = Fair - Fsg;
dWc    = (Fsc*Csc-FH) - (Fsp*Crgc+12.0*Fair*(cocyc+co2cyc));
dTreg  = (Qin - Qout)/((Wreg + Wsp)*cpc + MI); 

%>>>>>>>Reactor Riser Energy Calculations'

dHcrak = 172.7  + 3.0*(Tr - Tref);
Qrcat  = Frgc*cpc*(Tr-Tbase);
Qslury = F4*(cpsv*(Tr-Tref)+Qsr);
Qff    = F3*(cpfv*(Tr-Tref)+Qfr);
Qcrak  = (F3 + F4)*dHcrak;
Qrin   = Qrgc + F3*cpfl*(T2-Tbasef);
Qrout  = Qrcat + Qslury + Qcrak + Qff;


dWr    = Frgc - Fsc;
dTr    = (Qrin - Qrout)/(Mcpeff);
dCsc   = (Frgc*Crgc + Fcoke - Fsc*Csc - Csc*dWr)/Wr;
dP5    = 0.833*(Fwg - FV11 - FV12 + FV13);

%>>>>>>>Power and electric motor AMP Calculations'
FlowCAB=(((F6/18.01528)*(10.731557)*(Tatm+459.67))/P1)*60;%ft^3/min
FlowWGC=(((Fwg/453.59)*(10.731557)*(560.96))/P7)*60;%ft^3/min
PCAB=(((Tatm+459.67)/0.9)*(14.7/491.67)*FlowCAB*12^2*((P2/P1) - 1))/33000;%HP
PWGC=((1/2)*(((Tcondenser-273.15)*1.8+491.67)/0.9)*(14.7/491.67)*FlowWGC*12^2*(Crw^2 - 1.0d0))/33000;%HP
ACAB=PCAB*746/(13800*0.9*1.1*1.73);%Amp
AWGC=PWGC*746/(13200*0.9*0.98*1.73);%Amp

%>>>>>>>derivatives'
delta(1) = (F5*DHfu-UAf*Tlm-Qloss)/taufb;
delta(2) = (T2ss-T2)/taufo;
delta(3) = 5.0*(FV11 - F11);
delta(4) = dP5;
delta(5) = (R*(Tdla+Tf)/29.0/Vdla)*(F8-FV8-F9-F10);
delta(6) = R*(Rn*dTreg+(Treg+Tf)*dn)/Vregg;
delta(7) = ((Fsc/vcatl/Alp)+rhoag-rho)/taufil;
delta(8) = R*(Tdca+Tf)*(F6-FV7-F7)/29.0/Vdca;
delta(9) = dCsc;
delta(10)= (dWc - Crgc*(Fsc-Fsp))/Wreg;
delta(11)= dTreg;
delta(12)= Fsp - Frgc;
delta(13)= Fsc - Fsp;
delta(14)= dn;
delta(15)= dWr;
delta(16)= dTr ;      
delta(17)= (-Fair+(F7+F9+F10)/29.0);
delta(18)= -Pblp+P6+(rho*hlift/144.0)+(zbed-zlp)*rhocdn/144.0;
delta(19)= R*(Tatm+Tf)*(Fv6-F6)/29.0/Vcms;
delta(20)=dPreHeatE;
delta(21)=dMainFracE;
delta(22)=dRegenE;
delta(23)=dRegenTE;
delta(24)=dCatIE;
delta(25)=dReacTE;
delta(26)=dGFpvgo;
delta(27)=dGFp1;
delta(28)=dGFp2;
delta(29)=dGFp3;
delta(30)=dGFp4;
delta(31)=dGFc5;
delta(32)=dGFb;
delta(33)=dGFp;
delta(34)=dGFe;
delta(35)=dGFm;
delta(36)=((Fwgi-Fwg))/5;
delta(37)=(Fcokei-Fcoke)/5;

%>>>>>>Define the output'

fac=xco*28+xco2*44+xo2*32+22.12;
yp(1)=P4;
yp(2)=deltP;
yp(3)=Fair*29.0d0;
yp(4)=P6;
yp(5)=T3;
yp(6)=T2;
yp(7)=Tr;
yp(8)=Treg;
yp(9)=splev;
yp(10)=Tcyc;
yp(11)=Tcyc-Treg;
yp(12)=xco*28*1.e+6/fac;
yp(13)=xo2*Fair*100/Fsg;
yp(14)=Csc;
yp(15)=Crgc;
yp(16)=Fair;
yp(17)=Wris;
yp(18)=Wreg;
yp(19)=Wsp;
yp(20)=P5;
yp(21)=V4;
yp(22)=V6;
yp(23)=V7;
yp(24)=V3;
yp(25)=V1;
yp(26)=V2;
yp(27)=T2;
yp(28)=P5;
yp(29)=P6;
yp(30)=Tr;
yp(31)=Treg;
yp(32)=Wr;
yp(33)=Frgc*60;
yp(34)=Fsp*60;
yp(35)=Fcoke*60;
yp(36)=yfinal(1);
yp(37)=yfinal(2);
yp(38)=yfinal(3);
yp(39)=yfinal(4);
yp(40)=yfinal(5);
yp(41)=yfinal(6);
yp(42)=yfinal(7);
yp(43)=yfinal(8);
yp(44)=yfinal(9);
yp(45)=yfinal(10);
yp(46)=ones(1,11)*yfinal-F3;%Reactor MB
yp(47)=ACAB;%Power CAB
yp(48)=AWGC;%Power WGC
yp(49)=F5*60;
yp(50)=F7*60;
yp(51)=Fsg*60;
yp(52)=FV11*60;

end