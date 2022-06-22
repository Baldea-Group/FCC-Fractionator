function [y,Zl,Keq,HL,HV]=Enthalpy(T,P,x,Mw)
%Computes the Enthalpies of the stage 
%Inputs:
%T-->Temperature of stage
%P-->Pressure of stage
%x-->Liquid composition of stage
%Mw-->Molecular weight of pseudo-components
%Outputs:
%y-->Vapor composition of stage
%Zl-->Liquid phase compressibility factor
%Keq-->Equilibrium constant
%HL-->Enthalpy liquid
%HV-->Enthalpy vapor

k=10;%# of components
fixc=4;%# of fixed components (Methane, Ethane, Propane, Butane)

%Computation of critical temperature, pressure and accentric fractor of pseudo-components
wpse=[];
Tcpse=[];
Pcpse=[];
for i=1:k-fixc
w=-5.55305029654881e-07*Mw(i)^2+0.00236172564924776*Mw(i)+0.0885182765841059;
T1=1.76502197695905e-11*Mw(i)^5-4.36256156911342e-08*Mw(i)^4+4.13720500729747e-05*Mw(i)^3-0.0193203380223777*Mw(i)^2+5.15466244955443*Mw(i)+185.407202849323;
P1=-3.94956174670886e-13*Mw(i)^5+1.18757248039115e-09*Mw(i)^4-1.40250698603397e-06*Mw(i)^3+0.000827965933144259*Mw(i)^2-0.264543405611834*Mw(i)+48.4998241245831;
wpse=[wpse;w];
Tcpse=[Tcpse;T1];
Pcpse=[Pcpse;P1];
end

%Methane
Tf1=190.6990051;%(K)
Pf1=46.40670577;%(bar)
wf1=0.01150;

%Ethane
Tf2=305.428009;%(K)
Pf2=48.83839996;%(bar)
wf2=0.09860;

%Propane
Tf3=369.8980103;%(K)
Pf3=42.56651352;%(bar)
wf3=0.15240;

%Butane
Tf4=425.1990051;%(K)
Pf4=37.96612264;%(bar)
wf4=0.20100;

%EOS parameters
R=8.31445/(100);%(barL/Kmol)
TcT=[Tf1;Tf2;Tf3;Tf4;Tcpse];
wT=[wf1;wf2;wf3;wf4;wpse];
PcT=[Pf1;Pf2;Pf3;Pf4;Pcpse];
 
ai=[];
bi=[];
b=0; 
as=0;
for i=1:k
Tr(i)=T/TcT(i);%Reduced Temperature
m(i)=0.480+(1.574*wT(i))-0.176*(wT(i)^2);
a1=(1+m(i)*(1-sqrt(Tr(i))))^2;%alpha(i)
A1=0.42747*(R*TcT(i))^2*a1/PcT(i);
ai=[ai;A1];
B1=0.08664*R*TcT(i)/PcT(i);
bi=[bi;B1];
b=b+x(i)*bi(i);%bm
as=as+x(i)*ai(i)^(0.5);
end
B=b*P/(R*T);

a=0;
for i=1:k
    for j=1:k
    Aij2=x(i)*x(j)*sqrt(ai(i)*ai(j));%(alpha)(am)
    a=a+Aij2;
    end
end
A=a*P/(R*T)^2;

%EOS solution, polinomial roots
poli=[(1),(-1),(A-B-1*(B^2)),(-A*B)];
r=real(roots(poli));

Zv=max(r);
Zl=min(r);

%Computation of gas composition
Omegav=[];%Fugacity vapor
Omegal=[];%Fugacity liquid
Keq=[];%Equilibrium constant
y=[];%Gas composition

for i=1:k
Omegav1=exp(((Zv-1)*(bi(i)/b))-log(Zv-B)-(a/(b*R*T))*((1/a)*(2*ai(i)^(0.5)*as)-(bi(i)/b))*(log((Zv+B)/(Zv))));
Omegav=[Omegav;Omegav1];
Omegav2=exp(((Zl-1)*(bi(i)/b))-log(Zl-B)-(a/(b*R*T))*((1/a)*(2*ai(i)^(0.5)*as)-(bi(i)/b))*(log((Zl+B)/(Zl))));
Omegal=[Omegal;Omegav2];
Keq=[Keq;(Omegal(i)/Omegav(i))];
y=[y;(x(i)*Keq(i))];
end

%Enthalpy of ig for pseudo-components
Hig=[];
for i=1:k-fixc
Hform=3.18673341454555e-12*Mw(i)^5-7.66793947146919e-09*Mw(i)^4+6.99563705038894e-06*Mw(i)^3-0.00299458018039734*Mw(i)^2-1.23174567119212*Mw(i)-43.2101898494864; %(KJ/mol)    
Acp=0.0406170356952307*Mw(i)-0.496280964146270;
Bcp=(0.596633232248030*Mw(i)+2.34536350241800)*(10^-3);
Ccp=(-0.191164134529219*Mw(i)-0.336824623327747)*(10^-6);
Hg1=((Acp*(T))+(Bcp*(T^2)/2)+(Ccp*(T^3)/3)-(Acp*(298.15))-(Bcp*(298.15^2)/2)-(Ccp*(298.15^3)/3))*(8.31445/(1000)); %KJ/mol
Hig1=Hform+Hg1;
Hig=[Hig,Hig1];
end

%Methane-Butane Enthalpy of formation
H1formfix= -74.9; %Methane (KJ/mol) 
H2formfix= -84.738; %Ethane (KJ/mol) 
H3formfix= -103.89; %Propane (KJ/mol) 
H4formfix= -126.19; %Butane (KJ/mol)

%Methane-Butane Cp(T) polynomial coefficients
A1cpf=1.702;%Methane
B1cpf=9.081*(10^-3);
C1cpf=-2.164*(10^-6);

A2cpf=1.131;
B2cpf=19.225*(10^-3);%Ethane
C2cpf=-5.561*(10^-6);

A3cpf=1.213;%Propane
B3cpf=28.785*(10^-3);
C3cpf=-8.824*(10^-6);

A4cpf=1.935;%Butane
B4cpf=36.915*(10^-3);
C4cpf=-11.402*(10^-6);

%Methane-Butane Cp(T) closed form integral 
Hg1f=((A1cpf*(T))+(B1cpf*(T^2)/2)+(C1cpf*(T^3)/3)-(A1cpf*(298.15))-(B1cpf*(298.15^2)/2)-(C1cpf*(298.15^3)/3))*(8.31445/(1000)); %Methane KJ/mol
Hg2f=((A2cpf*(T))+(B2cpf*(T^2)/2)+(C2cpf*(T^3)/3)-(A2cpf*(298.15))-(B2cpf*(298.15^2)/2)-(C2cpf*(298.15^3)/3))*(8.31445/(1000)); %Ethane KJ/mol
Hg3f=((A3cpf*(T))+(B3cpf*(T^2)/2)+(C3cpf*(T^3)/3)-(A3cpf*(298.15))-(B3cpf*(298.15^2)/2)-(C3cpf*(298.15^3)/3))*(8.31445/(1000)); %Propane KJ/mol
Hg4f=((A4cpf*(T))+(B4cpf*(T^2)/2)+(C4cpf*(T^3)/3)-(A4cpf*(298.15))-(B4cpf*(298.15^2)/2)-(C4cpf*(298.15^3)/3))*(8.31445/(1000)); %Butane KJ/mol

%Methane-Butane Hig
Hig1f=H1formfix+Hg1f; %Methane (KJ/mol)
Hig2f=H2formfix+Hg2f; %Ethane (KJ/mol)
Hig3f=H3formfix+Hg3f; %Propane (KJ/mol)
Hig4f=H4formfix+Hg4f; %Butane (KJ/mol)

HVTotal=[Hig1f,Hig2f,Hig3f,Hig4f,Hig]; %Vector of Hig

HLTotal=HVTotal;

HL1=0;
for i=1:k
  HL1=HL1+(HLTotal(i))*x(i)/1000; %Hl (KMJ/mol or 1000*KJ/mol) 
end

HV1=0;
for i=1:k
  HV1=HV1+(HVTotal(i))*y(i)/1000; %Hv (KMJ/mol or 1000*KJ/mol)
end

%Increment (forward finite difference approximation of df(T)/dT)
dT=0.01;
Td=dT+T;

aid=[];%Follows same definition as above but with Td
for i=1:k
Tr(i)=Td/TcT(i);
m(i)=0.480+(1.574*wT(i))-0.176*(wT(i)^2);
a1=(1+m(i)*(1-sqrt(Tr(i))))^2;
A1=0.42747*(R*TcT(i))^2*a1/PcT(i);
aid=[aid;A1];
end

ad=0;
for i=1:k
    for j=1:k
    Aij2=x(i)*x(j)*sqrt(aid(i)*aid(j)); 
    ad=ad+Aij2;
    end
end

HLa=(R*T*(Zl-1))+((1/b)*(log((Zl+B)/(Zl)))*(T*((ad-a)/dT)-a));%(Lbar/mol)
HVa=(R*T*(Zv-1))+((1/b)*(log((Zv+B)/(Zv)))*(T*((ad-a)/dT)-a)); %(Lbar/mol)

HL=HLa*(1/10000)+HL1;%(KMJ/mol or 1000*KJ/mol)
HV=HVa*(1/10000)+HV1;%(KMJ/mol or 1000*KJ/mol)


end