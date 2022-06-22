%Script that plots (the most relevant) results

%Results
out1=R7;                                                                                                                                                                                                                     
outss=SPfcc;
out2=R10;
out3=R6;
out4=R9;
SPT=SPfrac(:,2);
SPH1=SPfrac(:,1);
SPC1=SPfrac(:,3); %HN
SPC2=SPfrac(:,4); %LCO
out2ss=[SPH1,SPT,SPC1,SPC2];
out5=[Holdup(1,:)',(Temperature(2,:)'-273.15)*(9/5)+32,(Ttrack(1,:)'),(Ttrack(2,:)')];
out6=[ProductsMB(end,:)',ProductsMB(4,:)',ProductsMB(5,:)',ProductsMB(6,:)'];
out7=Valvesfrac;

%Axis labels
zfcc={'Tpre (F)';'Pfra (psig)';'Preg (psig)';'Trea (F)';'Treg (F)';'Lrea (lb)';};
MVfcc={'Ff (scf/min)';'Frgc (lb/min)';'Fsc (lb/min)';'Fa (lb/min)';'Ffg (mol/min)';'Flpg (mol/min)';};
Valvesfcc={'V4 (%)';'V6 (%)';'V7 (%)';'V3 (%)';'V1 (%)';'V2 (%)';};
VI={'CAB (Amp)';'WGC (Amp)';};
MVf={'Fr (lb/min)';'Fln (lb/min)';'Fhn (lb/min)';'Flco (lb/min)';};
Prod={'VGO (lb/min)';'LPG (lb/min)';'LN (lb/min)';'HN (lb/min)';'LCO (lb/min)';'Slurry (lb/min)'};
Valvesf={'V9 (%)';'V8 (%)';'V10 (%)';'V11 (%)'};
zf={'Lfra (kmol)';'Tfra (F)';'Thnt (F)';'Tlcot (F)';};
xx={'min'};

figure(1)
for i=1:6
loc=300+20+i;
subplot(loc)
if i==6
plot(out1(:,1),out1(:,i+1)*(1/1),out1(:,1),outss(:,i)*(1/1),'--k');ylabel(zfcc(i));xlabel(xx) 
else
plot(out1(:,1),out1(:,i+1),out1(:,1),outss(:,i),'--k');ylabel(zfcc(i));xlabel(xx)    
end

end

figure(2)
for i=1:6
loc=300+20+i;
subplot(loc)
plot(out2(:,1),out2(:,i+1));ylabel(MVfcc(i));xlabel(xx)    
end

figure(3)
for i=1:6
loc=300+20+i;
subplot(loc)
plot(out3(:,1),out3(:,i+1));ylabel(Valvesfcc(i));xlabel(xx)    
end

figure(4)
for i=1:2
loc=300+20+i;
subplot(loc)
plot(out4(:,1),out4(:,end-(2-i)));ylabel(VI(i));xlabel(xx)

end

figure(5)
for i=1:4
loc=300+20+i;
subplot(loc)
plot(time,out5(:,i),time,out2ss(:,i),'--k');ylabel(zf(i));xlabel(xx) 
end

figure(6)
for i=1:4
loc=300+20+i;
subplot(loc)
plot(time,out6(:,i));ylabel(MVf(i));xlabel(xx) 
end

figure(7)
for i=1:4
loc=300+20+i;
subplot(loc)
plot(time,out7(:,i));ylabel(Valvesf(i));xlabel(xx) 
end

figure(8)
for i=1:6
loc=300+20+i;
subplot(loc)
if i==1
plot(time,ProductsMB(i,:)');ylabel(Prod(i));xlabel(xx)    
else    
plot(time,ProductsMB(1+i,:)');ylabel(Prod(i));xlabel(xx) 
end
end
