function y=Filter(t,x,flag,u)
%Filter to connect FCC and Fractionator
%Inputs:
%x-->states (Filter)
%u-->FCC reactor outputs
%Outputs:
%y-->Filtered values

y=model(x,u);


end

function y=model(x,u)
%model-->First order filter 

y=[];
for i=1:13
 y1=(u(i)-x(i))/1800; 
 y=[y;y1];
end


end
