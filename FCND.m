function F=FCND(t,xfcc,flag,dist,ufcc,Flpg,Tcondenser,i)   
%Takes FCC function to run the ODE subroutine
%Inputs:
%xfcc-->Vector of states
%dist-->Vector of disturbances
%ufcc-->Vector of manipulated variables
%Flpg-->Flow of LPG
%Tcondenser-->Condenser temperature
%i-->Time step
%Outputs:
%F-->Derivatives
%yp-->FCC process "measurements" (global variable)

global yp

[F,yp]=FCC(xfcc,dist,ufcc,Flpg,Tcondenser,i);

F=F';

