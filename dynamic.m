function [R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,SPfcc,Temperature,Vapor,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Liquid,Holdup,EnthalL,EnthalV,Ttrack,SPfrac,ProductsMB,LPG,MVfrac,ufra,xfra,xc,products,errord,xfil,Valvesfrac,time]=dynamic()
%***********************************************************************************************************************************************************************************************************************************************************
%Script that runs the FCC-Fractionator (no inputs required)
%Different disturbances, set point values, MVs (e.g. pump around duties, process feed) and initial conditions can be established 
%Table 2 in manuscript provides the operating region and set point changes magnitudes. To avoid numerical problems (e.g. dry up column), it is recommended to wait a period of time (e.g. at least 30 minutes of simulation) after consecutive set point changes 
%(especially after applying the absolute maximum set point change for all controllers simultaneously. In this case, the waiting time for another set of set point changes should be more than 30 minutes of simulation)  
%If the set point for the controllers must be changed at each minute of simulation, the absolute maximum set point change value should be reduced (particularly for TC2,TC3,PC1,PC2). A "reduction" factor in the range [1/10,1/4] is recommended
%The intuition behind this scaling is that the process might traverse the operating region faster (since at each minute new set points are established) leading to numerical problems, so the moves are "suppressed" to avoid this issue
%The operating region and set point magnitudes described above are general guidelines. The user is welcome to try different set point values and magnitudes just be careful with potential numerical difficulties 
%ST-->SETS THE SIMULATION TIME
%Outputs:
%R1-R10--> "Measurements" from FCC. FCC.m defines the elements of each matrix
%SPfcc--> Set points FCC
%Temperature--> Temperature of fractionator stages (from top to bottom)
%Vapor--> Vapor flowrate of fractionator stages (from top to bottom)
%Comp1-Comp10--> Liquid composition of component(i, where i has 10 elements) for all stages (from top to bottom)
%Liquid--> Liquid flowrate of fractionator stages (from top to bottom)
%Holdup--> Holdup of fractionator stages (from top to bottom)
%EnthalL--> Enthalpy of liquid phase of fractionator stages (from top to bottom)
%EnthalV--> Enthalpy of vapor phase of fractionator stages (from top to bottom)
%Ttrack--> Cut points for products [HN,LCO]
%SPfrac--> Set points fractionator
%ProductsMB--> Feed, coke, products, slurry and reflux flowrates
%LPG--> LPG molar flowrate
%MVfrac--> MVs fractionator (duties, PA flowrates)
%ufra-->Decision variables for NLP (The vector is concatenated as a column: [Vapor flowrate; Temperature])
%xfra-->States (The vector is concatenated as a column: [Hold-up; Enthalpy liquid; Enthalpy vapor; Liquid flowrate])
%xc-->Matrix of liquid composition (Each row (i) of the matrix represents the liquid composition of the stage. Each column (j) represents the component composition)
%products-->Reflux,HN and LCO flowrates
%errord-->Integral error for fractionator controllers
%xfil-->Filtered reactor output values
%Valvesfrac-->Fractionator valve position value
%time-->time step
%***********************************************************************************************************************************************************************************************************************************************************
global yp
%%%%%%%%%%%%%%%%%%%%%%%%%%FCC Initial Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
xfcc=[1564.06758445685,616,24.0440758965397,24.9000000145994,40.0420669348231,28.0000000041903,3.16075658859897,35.0875282505755,0.0101439627122109,0.00250580481493606,1249.99997383664,9350.79072581430,271210.435908199,246.810877809621,98098.9991960623,968.999859177516,2.68045714195376,30.4463683645626,14.6380093682365,-24.2765676899209,-0.612986761842691,0.611929610398922,-5641.47052429781,-30153150.1669079,-1879.52204006479,0.00468987379029097,0.0481251434795964,0.527273962836787,0.357213788281919,2.47339657486076,1.41455228053760,2.70179775410060,2.44381704476869,0.631637205219422,0.627951527686494,457.718461235902,6.31203102548586];
ufcc=[165,0,1,0,0,0,0.431260000000000,0,0,616,24.9000000000000,28,1250,98100,969];% LAST 6 COMPONENTS ARE THE SET POINTS FOR FCC. 
%Detailed definiton of each ufcc component can be found in FCC.m file
dist=[75,25,460.900000000000,0.900000000000000]; %DISTURBANCES
%Detailed definiton of each dist component can be found in FCC.m file
Flpg=457.718467816417;
Tcondenser=310.707095521277;

%%%%%%%%%%%%%%%%%%%%%Fractionator Initial Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MV=[0.747045932872508;298.895742822403;32;216.5;35.0441537396814;1243.39471672084];
ufra=[1647.78648413911;3547.63742893480;3382.54085689940;3077.23029277231;3039.30039378578;3038.30153796675;3031.04170871230;2992.96394430377;2929.38471956076;3035.21348344919;3023.22094129519;3055.92931441955;3054.01364490236;3044.12036102462;3036.53284553629;3032.88546649540;3031.43608755138;3030.73615977709;3031.85542990300;3317.69553572937;310.707095521277;392.000195034657;423.489189888333;449.071037007143;455.829436826681;457.278865297285;458.591772080347;463.785597803382;472.370698542057;509.816789386916;517.590163820810;519.337427827793;520.027072886078;521.878006501369;523.251988845496;523.948479342715;524.270766239095;524.490435402004;594.289627412580;628.139421202578];
xfra=[69.9992488400148;10.6316859735314;5.54317455602710;4.91100895722520;4.89436078215364;3.09047029825112;2.45584074341618;1.39618700941233;3.74406892836788;3.54419293748662;3.50526289424775;3.47333455975866;0.581208632217122;0.454748252715454;0.393957106911256;0.369799493313600;0.358132864975684;21.1000283535949;25.8640187364202;0.332060339295454;-0.224208017297567;-0.246522780582604;-0.310253017194136;-0.335719668638253;-0.338098963445289;-0.338858427183761;-0.345902800895044;-0.377049332312538;-0.414759572319589;-0.430408083956499;-0.433430113589529;-0.434099943180356;-0.435987977736354;-0.447571182665404;-0.455658696962468;-0.459238836223210;-0.460571288645674;-0.461449885379717;-0.438370329521518;-0.445316469740412;-0.112072792215712;-0.142025512847148;-0.142222571040525;-0.138762351945039;-0.139091175556689;-0.139256499719380;-0.138852001804254;-0.136654245646270;-0.133003418846954;-0.136537633847564;-0.136534687842125;-0.136518050403360;-0.136377516308497;-0.135533047296566;-0.134886368848381;-0.134573362425308;-0.134446050237712;-0.134380488991606;-0.174913487206278;-0.188804971011267;802.998070714435;637.901158411887;332.590473361626;294.660537433512;293.661646929218;185.428217895067;147.350444604971;83.7712205647395;224.644135702073;212.651576249197;210.315773654865;208.400073585520;34.8725179330273;27.2848951629272;23.6374264146753;22.1879695988160;21.4879718985410;1266.00170121569;1551.84112418521;16.6030169647727];
xc=[0.000730235602087646,0.00372006166875014,0.0449160106227330,0.139356343425088,0.240492627626492,0.569407502810893,0.00137721823703024,6.92087648762089e-12,2.04168509426954e-26,-2.61959392296970e-32;0.000334414014146372,0.00101897523290478,0.00909809956114959,0.0245075757206946,0.0715197965665978,0.805646930325290,0.0878741953189763,1.32602379525904e-08,1.29974453508235e-21,-1.98249322448692e-26;0.000342118151941009,0.000893370382498736,0.00695902081404365,0.0156121143250025,0.0314699456549243,0.399614511665990,0.545107423462561,1.49554304830301e-06,2.96879527394061e-18,-2.75777850354419e-22;0.000371388883508717,0.000877397119630828,0.00634125405177895,0.0131255127127142,0.0215961342254390,0.171830359349128,0.785828872994828,2.90806629994151e-05,8.42206265615051e-16,-1.80720718410264e-19;0.000376010088946707,0.000867105833276753,0.00615685670928431,0.0125119274032362,0.0197775706298387,0.137004471256256,0.822924231092917,0.000381826986099046,1.49185925868036e-13,-5.52729229133610e-18;0.000377686686936986,0.000866438571228639,0.00612870663488887,0.0124061114677264,0.0194546642445798,0.131871471895551,0.824176152207220,0.00471876826742203,2.44541519224196e-11,-5.00042931705508e-18;0.000382424464862341,0.000872714357434769,0.00614904136734986,0.0123981956806917,0.0192966364602325,0.129311838919261,0.777598256545508,0.0539908885592787,3.64539017432342e-09,2.64911772315457e-18;0.000400024361159347,0.000894464046521220,0.00620717285729142,0.0123246547582017,0.0186353777333138,0.120111722966723,0.567507352385685,0.273919017903400,2.12987711525807e-07,3.97613996517875e-15;0.000427241272341516,0.000925755578058187,0.00627589413791006,0.0121702193459629,0.0176002587126851,0.106771477968750,0.303198883002485,0.552625672209204,4.59777197822388e-06,6.31957010925434e-13;0.000407916325876747,0.000793998858013799,0.00498411857163051,0.00895451675892823,0.0112995627973560,0.0588093259148082,0.131227318685345,0.783477776273698,4.54657911118516e-05,2.32345765950015e-11;0.000412498322946457,0.000786765714755351,0.00486595642129017,0.00861093864290391,0.0105525371159684,0.0525248500042417,0.0817559529664807,0.840138634350095,0.000351865520279224,9.41039657478064e-10;0.000415155511979851,0.000788232479792742,0.00485884617664370,0.00856952430654025,0.0104356385084450,0.0514852009202325,0.0715441547694136,0.849372900295214,0.00253031204434431,3.49873936726212e-08;0.000418360041155387,0.000792769986141187,0.00487966144310618,0.00859353569841551,0.0104360573494958,0.0513038373319371,0.0692294158287375,0.836734392764956,0.0176107175926509,1.25196340150351e-06;0.000427546260085840,0.000805887972276076,0.00494005727792305,0.00866390904419732,0.0104402136533784,0.0508309412334667,0.0664551855722971,0.748314031869416,0.109083000839452,3.92262775067903e-05;0.000435251299540706,0.000817153698635069,0.00499374213634322,0.00873095847947023,0.0104602272900075,0.0505628114883548,0.0647351428226610,0.685024872680089,0.174050512715489,0.000189327389406824;0.000440155250103753,0.000824648894646010,0.00503150547762483,0.00878283682022606,0.0104908347464279,0.0505221650529794,0.0640013231762480,0.656191870114909,0.203040134656338,0.000674525810499191;0.000443502100289954,0.000830064244467471,0.00506053522234762,0.00882646341733982,0.0105273010849871,0.0506049588009945,0.0637876034966238,0.645424036922622,0.212336574881648,0.00215895982867887;0.000446517466546717,0.000835077087487109,0.00508812112220212,0.00886935714915974,0.0105668971603744,0.0507273893435450,0.0637247509928538,0.640446001911348,0.212719562435708,0.00657632533077642;0.000332071030021883,0.000539365817469141,0.00297548888672389,0.00470816869014084,0.00476067446478888,0.0195481864673815,0.0266278610329684,0.544238039047297,0.379906302722132,0.0163638418410765;0.000320842222488389,0.000491881890237686,0.00259910319544938,0.00393312933850515,0.00364120715415336,0.0130801198927189,0.0112100934166477,0.281216487377775,0.607390131596276,0.0761170039157492];
SP=[70;245.93;530.33;755.33]; %SET POINTS FOR FRACTIONATOR
Distillateini=xfra(20*3+1)/MV(1); 
products=[Distillateini;100.973572040226;163.634173822281];
errord=[0.136692017830002;-0.0707331131770190;-0.352870002514544;-2.77937044322705];
Xfilin=[3025.85220559869,793.705291729987,0.0559148422896900,0.0562430796039667,0.217605876121722,0.240577895792075,0.125956811268771,0.220240485954903,0.0318076418832473,0.0469505105292196,0.00428524729043614,0.000417609265968941,1.71679455003844]'; %States for the filter (dist+P)

%%%%%%%%%%%%%%%%%%%%%%FCC VARIABLES(RESULTS) DEFINITION%%%%%%%%%%%%%%%%%%%%
[delta,yp]=FCC(xfcc,dist,ufcc,Flpg,Tcondenser,1);

T=0.0;
R1=[T  yp(1) yp(2) yp(3) yp(4) yp(5)];
R2=[T  yp(6) yp(7) yp(8) yp(9) yp(10)];
R3=[T  yp(11) yp(12) yp(13) yp(14) yp(15)];
R4=[T  yp(16) yp(17) yp(18) yp(19) yp(20)];
R5=[T  xfcc];
R6=[T  yp(21) yp(22) yp(23) yp(24) yp(25) yp(26)];
R7=[T  yp(27) yp(28) yp(29) yp(30) yp(31) yp(32)];
R8=[T  yp(33) yp(34) yp(35)];
R9=[T  yp(36) yp(37) yp(38) yp(39) yp(40) yp(41) yp(42) yp(43) yp(44) yp(45) yp(46) yp(47) yp(48)]; 
R10=[T  yp(49) yp(33) yp(34) yp(50) yp(51) yp(52)]; 
SPfcc=[ufcc(10) ufcc(11) ufcc(12) ufcc(15) ufcc(13) ufcc(14)];
options=[];
h = 10; %%%LOOK FOR TIME STEP%%% (10 seconds)

%%%%%%%%%%%%%%%%%%%%FRACTIONATOR VARIABLES(RESULTS) DEFINITION%%%%%%%%%%%%%
Temperature=[];
Vapor=[];
Comp1=[];
Comp2=[];
Comp3=[];
Comp4=[];
Comp5=[];
Comp6=[];
Comp7=[];
Comp8=[];
Comp9=[];
Comp10=[];
Liquid=[];
Holdup=[];
EnthalL=[];
EnthalV=[];
Ttrack=[];
time=[];
SPfrac=[];
Valvesfrac=[];
ProductsMB=[];
LPG=[];
xfil=[];
MVfrac=[]; 
ST=500;%%%%%%%%%%%%%%%SET THE SIMULATION TIME(ST) IN MINUTES%%%%%%%%%%%%%%
for minute=1:ST%Iterate FCC-Fractionator solution over the simulation time%
    minute
      
% Set point tracking case study
 if (minute>=60)&&(minute<240)
     ufcc(10)=618.5;
 elseif (minute>=240)
     ufcc(15)=971;  
 end

%Disturbance rejection case study
% if (minute>=60)
%     dist(2)=20;
% end
  
for j=1:(60/h)%%%%%%%%%%%%%% INTEGRATE FCC-Fractionator %%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%FCC SOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tspan=[T T+h];
	[t,y]=ode15s('FCND',tspan,xfcc,options,dist,ufcc,Flpg,Tcondenser,minute);	 
	
      
   [n m]=size(y);
    xfcc=y(n,:);
    
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

    Frout=[(yp(36)/MWp5) (yp(37)/MWp4) (yp(38)/MWp3) (yp(39)/MWp2) (yp(40)/MWp1) (yp(41)/MWc5) (yp(42)/MWb) (yp(43)/MWp) (yp(44)/MWe) (yp(45)/MWm)]*(453.59); %(mole/s)
    FRT=Frout*ones(10,1);%(mole/s)
    xfeedfrac=(fliplr(Frout))/FRT;
    FtotalF=FRT*(60*60/1000); %Kmol/h
    Pin=yp(28)*(1/14.503773773);%bar
    ToutF=(yp(7)-32)*(5/9)+273.15;%K
    Dist=[FtotalF;ToutF;xfeedfrac'];
    
    ufilter=[Dist;Pin];
    [t1,yfilter]=ode15s('Filter',tspan,Xfilin,options,ufilter);
    [n m]=size(yfilter);
    Xfilin=yfilter(n,:); %Filtered data to Fractionator
    T=T+h;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%FRACTIONATOR SOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iswitch=rem(j,2);
    
if iswitch==1 %Reversed the solution order at each time step        
[Temperatureout,Vaporout,xcout,Liqout,Holdout,EnthalLout,EnthalVout,LN,HN,LCO,ELC2,ETC4,ETC5,ETC6,Ttrack1,Ttrack2,yout,Valvesf]=Fractionator(xfra,ufra,xc,MV,SP,products,errord,Xfilin,dist,minute);
else
[Temperatureout,Vaporout,xcout,Liqout,Holdout,EnthalLout,EnthalVout,LN,HN,LCO,ELC2,ETC4,ETC5,ETC6,Ttrack1,Ttrack2,yout,Valvesf]=Fractionatori(xfra,ufra,xc,MV,SP,products,errord,yout,Xfilin,dist,minute);    
end

xfra=[Holdout;EnthalLout;EnthalVout;Liqout];
ufra=[Vaporout;Temperatureout];
xc=xcout;
products=[LN;HN;LCO];
errord=[ELC2;ETC4;ETC5;ETC6];
MVfracout=[MV(1),MV(2),MV(3),MV(4),MV(5),MV(6)];
Flpg=(ones(1,10)*(yout(1,:)'*Vaporout(1)))*(1000/3600);%mol/s
Tcondenser=Temperatureout(1);
MWT=[MWm,MWe,MWp,MWb,MWc5,MWp1,MWp2,MWp3,MWp4,MWp5];
LPGMB=((1/3.6)*(1/453.59)*Vaporout(1)*yout(1,:).*MWT)*ones(10,1);
LNMB=((1/3.6)*(1/453.59)*LN*xcout(1,:).*MWT)*ones(10,1);
HNMB=((1/3.6)*(1/453.59)*HN*xcout(6,:).*MWT)*ones(10,1);
LCOMB=((1/3.6)*(1/453.59)*LCO*xcout(13,:).*MWT)*ones(10,1);
SMB=((1/3.6)*(1/453.59)*Liqout(end)*xcout(20,:).*MWT)*ones(10,1);
TMB=(ufcc(1)-(yp(35)/60+LPGMB+LNMB+HNMB+LCOMB+SMB))*(100/ufcc(1)); %Mass balance (%)
Conversion=((LPGMB+LNMB+HNMB+LCOMB)/ufcc(1))*100;%Conversion into valuable products (slurry not considered)
Liq1=((1/3.6)*(1/453.59)*Liqout(1)*xcout(1,:).*MWT)*ones(10,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FRACTIONATOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temperature=[Temperature,Temperatureout];
Vapor=[Vapor,Vaporout];
Comp1=[Comp1,xcout(:,1)];
Comp2=[Comp2,xcout(:,2)];
Comp3=[Comp3,xcout(:,3)];
Comp4=[Comp4,xcout(:,4)];
Comp5=[Comp5,xcout(:,5)];
Comp6=[Comp6,xcout(:,6)];
Comp7=[Comp7,xcout(:,7)];
Comp8=[Comp8,xcout(:,8)];
Comp9=[Comp9,xcout(:,9)];
Comp10=[Comp10,xcout(:,10)];
Liquid=[Liquid,Liqout];
Holdup=[Holdup,Holdout];
EnthalL=[EnthalL,EnthalLout];
EnthalV=[EnthalV,EnthalVout];
Ttrack=[Ttrack,[(Ttrack1-273.15)*9/5+32;(Ttrack2-273.15)*9/5+32]];%K
SPfrac=[SPfrac;SP'];
time=[time;minute];
MVfrac=[MVfrac;MVfracout];
ProductsMB=[ProductsMB,[ufcc(1);yp(35)/60;LPGMB;LNMB;HNMB;LCOMB;SMB;Liq1]*(60)]; %Feed,Coke,Gas,LN,HN,Die,Slurry,L1
LPG=[LPG,Flpg];
xfil=[xfil;Xfilin];
Valvesfrac=[Valvesfrac;Valvesf];

%%%%%%%%%%%%%%DISPLAY RELEVANT VARIABLES IN COMMAND WINDOW%%%%%%%%%%%%%%%%
TVFCC=table(yp(21),yp(22),yp(23),yp(24),yp(25),yp(26),'VariableNames',{'V4','V6','V7','V3','V1','V2'},'RowName',{'Valve(%)'});
disp(TVFCC)
TVF=table(Valvesf(1),Valvesf(2),Valvesf(3),Valvesf(4),'VariableNames',{'V9','V8','V10','V11'},'RowName',{'Valve(%)'});
disp(TVF)
TPF=table(LPGMB*60,LNMB*60,HNMB*60,LCOMB*60,SMB*60,'VariableNames',{'LPG','LN','HN','LCO','Slurry'},'RowName',{'Flowrate(lb/min)'});
disp(TPF)
TCMB=table(Conversion,TMB,'VariableNames',{'Conversion','Mass_Balance'},'RowName',{'%'});
disp(TCMB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FCC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R1=[R1; T/60  yp(1) yp(2) yp(3) yp(4) yp(5)];
R2=[R2; T/60  yp(6) yp(7) yp(8) yp(9) yp(10)];
R3=[R3; T/60  yp(11) yp(12) yp(13) yp(14) yp(15)];
R4=[R4; T/60  yp(16) yp(17) yp(18) yp(19) yp(20)];
R5=[R5; T/60  xfcc];
R6=[R6; T/60  yp(21) yp(22) yp(23) yp(24) yp(25) yp(26)];
R7=[R7;T/60  yp(27) yp(28) yp(29) yp(30) yp(31) yp(32)];  
R8=[R8;T/60  yp(33) yp(34) yp(35)];
R9=[R9;T/60  yp(36) yp(37) yp(38) yp(39) yp(40) yp(41) yp(42) yp(43) yp(44) yp(45) yp(46) yp(47) yp(48)];
R10=[R10;T/60  yp(49) yp(33) yp(34) yp(50) yp(51) yp(52)];
SPfcc=[SPfcc; ufcc(10) ufcc(11) ufcc(12) ufcc(15) ufcc(13) ufcc(14)];

end
 Plotall
end