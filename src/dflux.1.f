! NEED THIS HEADER INFORMATION

      subroutine dflux(flux,sol,kstep,kinc,time,noel,npt,coords,
     &     jltyp,temp,press,loadtype,area,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iscale,mi,
     &     sti,xstateini,xstate,nstate_,dtime)

      character*8 lakonl
      character*20 loadtype
!
      integer kstep,kinc,noel,npt,jltyp,konl(20),ipompc(*),nstate_,i,
     &  nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),node,idof,id,iscale,mi(*)
!
      real*8 flux(2),time(2),coords(3),sol,temp,press,vold(0:mi(2),*),
     &  area,co(3,*),coefmpc(*),sti(6,mi(1),*),xstate(nstate_,mi(1),*),
     &  xstateini(nstate_,mi(1),*),dtime
!
      intent(in) sol,kstep,kinc,time,noel,npt,coords,
     &     jltyp,temp,press,loadtype,area,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,mi,sti,
     &     xstateini,xstate,nstate_,dtime
     
      intent(out) flux,iscale

! END OF HEADER INFORMATION
! CUSTOM SUBROUTINE HERE
      
      x0=0
	y0=0 
	z0=0
      
      vx = 0.00
      vy = 0.02
      vz = 0.00
      
      xarc=vx*time(1)+x0
	yarc=vy*time(1)+y0
	zarc=vz*time(1)+z0
	  
      Xf=coords(1)-xarc	! coordinate of position x 
	Yf=coords(2)-yarc	! coordinate of position y 
	Zf=coords(3)-zarc	! coordinate of position z
	
	X=coords(1)
      Y=coords(2)
      Z=coords(3)

      IF(Yf.LE.Y) THEN
            flux(1)=50000
      ELSE
            flux(1)=0
      END IF

      !X=coords(1)
      !Y=coords(2)
      !Z=coords(3)
      
      !IF(X.LE.X0 .AND. X.GE.X1 .AND. Z.LE.Z0 .AND. Z.GE.Z1) THEN
      !      flux(1)=50000
      !ELSE
      !      flux(1)=0
      !END IF
      
      
      
      RETURN
      END 
