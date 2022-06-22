function zposition
%Integration of energy and MB along z axis
global Fair FH Treg 
global Tcyc o2cyc cocyc co2cyc
global xo2 xco xco2
		
      tz0=0.0;
      tz=tz0;
      zstep=0.5;

      xz(1)=Treg;
      xz(2)=(0.21*Fair-0.25*FH)/Fair;
      xz(3)=0.0;
      xz(4)=0.0;

      for i=1:4
			xzdot(i)=0.0;
      end
      
      while abs(tz-49.5)>1.0e-4
                  
         for i=1:4
            xzold(i)=xz(i);
         end 

         delta1=zcyc(tz,xz,xzdot);

         for i=1:4
		  		xz(i)=xzold(i)+zstep*delta1(i);
         end

         tz2=tz+zstep;
         delta2=zcyc(tz2,xz,xzdot);

         for i=1:4
            xz(i)=xzold(i)+zstep*(delta1(i)+delta2(i))/2.0;
         end 


         if abs(tz-45.0)>1.0d-4
            Tcyc=xz(1);
			o2cyc=xz(2);
			cocyc=xz(3);
			co2cyc=xz(4);
         end
         
         tz=tz+zstep;

      end
      	
      xo2=xz(2);
      xco=xz(3);
      xco2=xz(4);

      return





