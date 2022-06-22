function delta=zcyc(zz,xp,xpdot)
%Temperature profile and MB
global Fair FH vs rhocdl Crgc zbed P6 epse 

Areg=590.0;
dh1=46368.0;
dh2=169080.0;
cpco2=11.0;
cpc=0.31;
cpco=7.28;
cph2o=8.62;
cpn2=7.22;
cpo2=7.62;

zcyc=45.0;
Me=Areg*vs*rhocdl;

if(zz<zcyc)
  delz=1.0d0;
else
  delz=0.0d0;
end

cpz= 0.79*cpn2+xp(3)*cpco+xp(4)*cpco2+xp(2)*cpo2+(0.5*cph2o*FH+delz*cpc*Me)/Fair;

k1=6.95470*exp(19.88-34000.0/(xp(1)+459.6));
k2=0.69148*exp(15.06-25000.0/(xp(1)+459.6));
k3=0.6412*P6*exp(25.55-45000.0/(xp(1)+459.6));

if(zz<=zbed); rhob=1.0-epse; end
if(zz>zbed);  rhob=(1.0-epse)*exp(-1000.0*Fair/Areg/vs/rhocdl*(zz-zbed)); end
if(zz>zcyc);  rhob=0.0d0; end


 dpx2=(100.0*(-0.5*k1-k2)*rhob*Crgc-k3*xp(3))*xp(2)/vs;
 dpx3=(100.0*k1*rhob*Crgc-2*k3*xp(3))*xp(2)/vs;
 dpx4=-(dpx2+0.5*dpx3);

 if(zz<=zbed)
     delta(1)=-xpdot(1);
 else
     delta(1)=-xpdot(1)+(dh1*dpx3+dh2*dpx4)/cpz;
 end

 delta(2)=-xpdot(2)+dpx2;
 delta(3)=-xpdot(3)+dpx3;
 delta(4)=-xpdot(4)+dpx4;

 return


