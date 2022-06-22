
function [y,beta]=PFR(a0,F3,F4,vris,disturbance)
global Tr P4
x0=0;%Initial
xf=1;%Final
h=0.1;%Step
n=(xf-x0)/h;% # of steps

for i=1:n %Runge Kutta (Fourth order)
    y1=riser(a0,vris,disturbance);
    K1=y1*h;
    y2=riser((a0+K1/2),vris,disturbance);
    K2=y2*h;
    y3=riser((a0+K2/2),vris,disturbance);
    K3=y3*h;
    y4=riser((a0+K3),vris,disturbance);
    K4=y4*h;
    a0=a0+(K1/6)+(K2/3)+(K3/3)+(K4/6);
    yt=a0;
end
y=yt; %(mol j/Kg gas)

T=(Tr-32)*(5/9)+273.15;%K
P=P4*(1/14.5038);%bar
Rgas=8.314/100000;%(m^3bar/molK)
Ftotal=(F3+F4)*(1/2.20462); %(Kg/s)
Astripper=60;%(ft^2);
Lstripper=15;%ft 
Vstripper=Astripper*Lstripper*((1/3.28084)^3);%m^3

Mwave=1/(ones(1,10)*yt(1:end-1)); %Average Mw

beta=(Rgas*T*Ftotal)/(P*Vstripper*Mwave); %time constant for CST
 end


function y=riser(a,vris,disturbance)
global Tr Wris Frgc P4
%a concentration of j (mol j/Kg gas)
APInominal=25;

%Parameters MW (g/mol) lump (j)
MWp5=446.094358949534;%(g/mol)
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
MWT=[MWp5;MWp4;MWp3;MWp2;MWp1;MWc5;MWb;MWp;MWe;MWm;MWc];

%Stoichiometric Paremeters
v=[];
for i=1:9
    v1=[];
    if i<=5
    for j=i+1:11
    a1=MWT(i)/MWT(j);
    v1=[v1;a1];
    end
    else
    for j=i+1:10
    a1=MWT(i)/MWT(j);
    v1=[v1;a1];   
    end    
    end
    v=[v;v1];
    clear v1
end
a1=length(v);

%Kinetic information
%Pre-exponential Factor (frequency factor) (m^3/(Kgcath)). 
A=[685381.702376125;16071.7417336359;55910.6460564503;5340.56954492411;4261075.43228439;3702886.86760602;4314348.37218425;4221817.31608134;2594270.14855151;3439792.76335786;4298171.84154773;4309482.15658185;4188204.43659805;4171507.14955002;4268930.21102278;4235978.25533608;4035248.96093414;4159694.61983738;4016315.59411039;4317309.36319068;3998905.55010847;3977645.15248325;3946590.54946382;4093321.37808001;3768025.72738946;3473846.21400886;3710947.31883704;4249525.78762700;4168288.98051169;4201041.80321129;4198160.18799487;4004640.42044027;3732579.74974075;3597877.86970591;15009.5521341976;6922.30720980437;19041.3126135312;488155.088910786;2083212.12299930;61865.7674159749;4295780.23091937;4264499.48678023;4117346.09044240;3854696.38239672;4254615.63253292;4099182.31561032;4074269.75010443;4034030.70015591;4125807.31083366;4151745.82149610];
%Activation Energy (KJ/mol)
R=8.314/1000;%KJ/molK (Gas Constant)
E=[82.8791270746213;45.6320226868409;61.1861084125379;36.6372841257022;88.1355579147187;83.4702620995338;93.4587040123872;110.518807939558;114.747751498962;93.1736077915559;100.154816005713;99.4119449292334;100.092458367997;99.4111112794267;99.6338215122641;100.140463287950;102.298311600346;101.750324286216;110.028728198801;98.3105495562886;99.6339220645355;97.3374067932750;110.484273942267;101.009077284125;108.879891812183;131.360134407370;121.413630063829;100.849697879039;100.675404066808;101.928965649227;101.270128345124;106.053632596030;108.460949734182;118.240920291922;68.1750414644057;80.6231361510246;69.3215293814854;109.546981921518;132.564043129576;105.830211294964;102.742794601642;101.318881654179;107.392472585985;110.240875298551;95.7276754466431;108.735606818799;111.689214350240;106.288490074524;120.231269015535;97.8072408305545];

%Rate Value
T=(Tr-32)*(5/9)+273.15;%K
k=[];
for i=1:a1
cte=A(i)*exp(-1*(E(i))/(R*T)); 
k=[k;cte];    
end

%Feed disturbance(API) effect (a subset of the reactions only)
%Gas production%
k(5)=k(5)*(disturbance/APInominal)^0.15;
k(6)=k(6)*(disturbance/APInominal)^0.15;
k(7)=k(7)*(disturbance/APInominal)^0.15;
k(8)=k(8)*(disturbance/APInominal)^0.15;
k(9)=k(9)*(disturbance/APInominal)^0.15;

%Coke production%
k(10)=k(10)*(disturbance/APInominal)^(-0.15);

K=[ -1*(k(1)+k(2)+k(3)+k(4)+k(5)+k(6)+k(7)+k(8)+k(9)+k(10)),0,0,0,0,0,0,0,0,0,0
    v(1)*k(1),-1*(k(11)+k(12)+k(13)+k(14)+k(15)+k(16)+k(17)+k(18)+k(19)),0,0,0,0,0,0,0,0,0
    v(2)*k(2),v(11)*k(11),-1*(k(20)+k(21)+k(22)+k(23)+k(24)+k(25)+k(26)+k(27)),0,0,0,0,0,0,0,0
    v(3)*k(3),v(12)*k(12),v(20)*k(20),-1*(k(28)+k(29)+k(30)+k(31)+k(32)+k(33)+k(34)),0,0,0,0,0,0,0
    v(4)*k(4),v(13)*k(13),v(21)*k(21),v(28)*k(28),-1*(k(35)+k(36)+k(37)+k(38)+k(39)+k(40)),0,0,0,0,0,0
    v(5)*k(5),v(14)*k(14),v(22)*k(22),v(29)*k(29),v(35)*k(35),-1*(k(41)+k(42)+k(43)+k(44)),0,0,0,0,0
    v(6)*k(6),v(15)*k(15),v(23)*k(23),v(30)*k(30),v(36)*k(36),v(41)*k(41),-1*(k(45)+k(46)+k(47)),0,0,0,0
    v(7)*k(7),v(16)*k(16),v(24)*k(24),v(31)*k(31),v(37)*k(37),v(42)*k(42),v(45)*k(45),-1*(k(48)+k(49)),0,0,0
    v(8)*k(8),v(17)*k(17),v(25)*k(25),v(32)*k(32),v(38)*k(38),v(43)*k(43),v(46)*k(46),v(48)*k(48),-1*(k(50)),0,0
    v(9)*k(9),v(18)*k(18),v(26)*k(26),v(33)*k(33),v(39)*k(39),v(44)*k(44),v(47)*k(47),v(49)*k(49),v(50)*k(50),0,0
    v(10)*k(10),v(19)*k(19),v(27)*k(27),v(34)*k(34),v(40)*k(40),0,0,0,0,0,0];
    
%Inert Aromatic Function
Kh=0.128;
Car=0;% aromatcis in feed (%wt)
fcar=1/(1+Kh*Car);

%Basic Nitrogen Adsorption
Kn=0.01; %adsorption coefficient (g nitrogen/g catalyst)
N=0;% grams of nitrogen in catalyst
gN=1/(1+Kn*N);

%Catalyst decay (Pressure dependance can be incorporated)
beta=470.259556216469; 
gama=0.652185139926354;
tc=(Wris/Frgc)*(1/3600); % residence time catalyst riser (hours)
phi=1/(1+(beta*(tc^gama)));

%%%Rate expression%%%
rsarea=9.6; %ft^2 (Aris in McFarlane study, cross sectional riser area)
vsup=vris/rsarea;%ft/s
epsilon=min(1,(0.332+0.06*vsup)); % Used void calculations employed in Regenerator


Rgas=8.314/100000;%(m^3bar/molK)
P=P4*(1/14.5038);%bar
densityg=P/(Rgas*T*(ones(1,11)*a));%(Kg/m3)
densityc=45*(3.28084^3/2.20462);%(Kg/m3) (Used density in catalyst coming from U-Bents)
Gv=((vris-(Frgc/68))*(3600/3.28084^3))*densityg/(rsarea*(1/3.28084^2));%Kg gas/m2h
Swh=(Gv*epsilon)/((60*(1/3.28084))*densityc);%(Kggas/Kgcat-h) hris=60 ft (heigth of riser)

y=(fcar)*(gN)*(phi)*(P/(Swh*Rgas*T))*(1/(ones(1,11)*a))*(K*a);

end
