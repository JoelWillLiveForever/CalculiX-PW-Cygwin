!     NEED THIS HEADER INFORMATION
!     Abaqus interface
!     subroutine dflux(flux,sol,kstep,kinc,time,noel,npt,coords,
!     1 jltyp,temp,press,sname)
      
!     Calculix interface 
      subroutine dflux(flux,sol,kstep,kinc,time,noel,npt,coords,
     &     jltyp,temp,press,loadtype,area,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iscale,mi,
     &     sti,xstateini,xstate,nstate_,dtime)

      implicit real*8(a-h,o-z)
      
      DIMENSION COORDS(3),FLUX(2),TIME(2)
      
      character*80 sname     
      
      character*8 lakonl
      character*20 loadtype

!     intent(in) sol,kstep,kinc,time,noel,npt,coords,
!     &     jltyp,temp,press,loadtype,area,vold,co,lakonl,konl,
!     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,mi,sti,
!     &     xstateini,xstate,nstate_,dtime
      
!     intent(out) flux,iscale

!     END OF HEADER INFORMATION
!     CUSTOM SUBROUTINE HERE
      
      pi=acos(-1.0)
      
      k=(kstep+2)/3
      
      call heat_line_input(k,nwps,jshape,rad,power,speed,eff,thold,
     &     a,b,cf,cr,x0,y0,z0,x1,y1,z1)

      timenet=time(1)-thold
      
      if(jshape.eq.1) then
         xarc=  x0
         yarc=  y0
         zarc= speed*timenet + z0
         xcod=coords(1)-xarc	
         ycod=coords(2)-yarc	
         zcod=coords(3)-zarc	      
      else if(jshape.eq.2) then
         omega=speed*time(1)/rad	
         ythe=cos(omega)*rad
         zthe=sin(omega)*rad	
         xthe=x0	
         xcod=coords(1)-xthe	
         ycod=coords(2)-ythe	 
         zcod=coords(3)-zthe	
      endif

!     write(7,*) "xcod,ycod,zcod",xcod,ycod,zcod

      call dellipse(power,eff,a,b,cf,cr,xcod,ycod,zcod,qflux)
      flux(1)=qflux    
      
      RETURN
      END
!------------------------------------------------------------------------------
!     Read dflux inputs
!     
      subroutine heat_line_input(k,nwps,jshape,rad,power,speed,eff,thold,
     &     a,b,cf,cr,x0,y0,z0,x1,y1,z1)

      implicit real*8(a-h,o-z)

      if(k.eq.1) then
         nwps=1
         jshape=2
         rad=64.894	
         power=1674.000000
         speed=1.693000
         eff=0.75    
         thold=0.25
         a=1.872644
         b=1.872644
         cf=1.872644
         cr=3.745288
         x0=0.00373 
         y0=0.0
         z0=64.894329
         x1=0.00373 
         y1=0.0
         z1=64.894329   
      elseif(k.eq.2) then
         nwps=1
         jshape=2
         rad=66.799	
         power=1710.000000
         speed=1.693000
         eff=0.80    
         thold=0.25
         a=6.099836
         b=0.9
         cf=6.099836
         cr=12.199671
         x0=0.008756 
         y0=0.0
         z0=66.798856
         x1=0.008756 
         y1=0.0
         z1=66.798856                        
      elseif(k.eq.3) then
         nwps=1
         jshape=2
         rad=69.958	
         power=1674.0
         speed=1.693000
         eff=1.05    
         thold=0.25
         a=5.770694
         b=5.770694
         cf=5.770694
         cr=11.541387
         x0=3.515449 
         y0=0.0
         z0=69.957824 
         x1=3.515449 
         y1=0.0
         z1=69.957824 
      elseif(k.eq.4) then
         nwps=1
         jshape=2
         rad=70.119	
         power=1674.0
         speed=1.693000
         eff=1.05    
         thold=0.25
         a=5.267341
         b=5.267341
         cf=5.267341
         cr=10.534682
         x0=-3.604096 
         y0=0.0
         z0=70.119194
         x1=-3.604096 
         y1=0.0
         z1=70.119194
      elseif(k.eq.5) then
         nwps=1
         jshape=2
         rad=74.036	
         power=1953.0
         speed=1.693000
         eff=0.84    
         thold=0.25
         a=4.788613
         b=4.788613
         cf=4.788613
         cr=9.577227
         x0=3.932621 
         y0=0.0
         z0=74.035517 
         x1=3.932621 
         y1=0.0
         z1=74.035517                    
      elseif(k.eq.6) then
         nwps=1
         jshape=2
         rad=73.879	
         power=2256.0
         speed=1.693000
         eff=0.65    
         thold=0.25
         a=4.562056
         b=4.562056
         cf=4.562056
         cr=9.124112
         x0=-4.428784 
         y0=0.0
         z0=73.878978
         x1=-4.428784 
         y1=0.0
         z1=73.878978 
      elseif(k.eq.7) then
         nwps=1
         jshape=2
         rad=76.35	
         power=2256.0
         speed=1.693000
         eff=1.33    
         thold=0.25
         a=7.817115
         b=1.25
         cf=7.817115
         cr=15.634231
         x0=-0.01223 
         y0=0.0
         z0=76.349887
         x1=-0.01223 
         y1=0.0
         z1=76.349887 
      endif 
      
      return      
      end      
!------------------------------------------------------------------------------
      subroutine dellipse(powerk,effk,ak,bk,cfk,crk,
     &     xcod,ycod,zcod,qflux)

      implicit real*8(a-h,o-z)
      
      pi=acos(-1.0)
      
      cf=cfk
      cr=crk
      if(zcod.ge.0.0) then
         ff=2*cf/(cf+cr)
         fheat=ff
         cheat=cf
      else if(zcod.lt.0.0) then
         fb=2*cr/(cf+cr)
         fheat=fb
         cheat=cr
      else
         stop "*** error in fheat calculation ***"
      end if

      qtop=6.0*sqrt(3.0)*effk*powerk*fheat      	
      qbot=pi*ak*bk*cheat*sqrt(pi) 
      qmax=(qtop/qbot) 

      gcutk=ak/10 
      if(abs(xcod).lt.gcutk) then
         xcod=gcutk
      end if
      if(abs(ycod).lt.gcutk) then
         ycod=gcutk
      endif
      if(abs(zcod).lt.gcutk) then
         zcod=gcutk
      endif        
      
      rx2=(xcod*xcod)/(ak**2)	
      ry2=(ycod*ycod)/(bk**2)
      rz2=(zcod*zcod)/(cheat**2)

      qflux=qmax*exp(-3.0*(rx2+ry2+rz2))
      
c     write(7,*) "zcod",zcod
c     write(7,*) "qtop,bbot,qmax=",qtop,qbot,qmax
c     write(7,*) "abc f ",ak,bk,cheat,fheat
c     write(7,*) "q =",qqflux

      return
      end
!------------------------------------------------------------------------------
