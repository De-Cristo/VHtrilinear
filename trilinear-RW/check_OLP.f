      PROGRAM DRIVER
C     *****************************************************************
C     ********
C     THIS IS THE DRIVER FOR CHECKING THE STANDALONE MATRIX ELEMENT.
C     IT USES A SIMPLE PHASE SPACE GENERATOR
C     *****************************************************************
C     ********
      IMPLICIT NONE
C     
C     CONSTANTS  
C     
      REAL*8 ZERO
      PARAMETER (ZERO=0D0)

      LOGICAL READPS
      PARAMETER (READPS = .FALSE.)

      INTEGER NPSPOINTS
      PARAMETER (NPSPOINTS = 4)

C     integer nexternal and number particles (incoming+outgoing) in
C      the me 
c      INTEGER NEXTERNAL, NINCOMING
c      PARAMETER (NEXTERNAL=7,NINCOMING=2)
      INCLUDE 'nexternal.inc'
      CHARACTER(512) MADLOOPRESOURCEPATH

C     
C     INCLUDE FILES
C     
C     the include file with the values of the parameters and masses   
C        
      INCLUDE 'coupl.inc'
C     particle masses
      REAL*8 PMASS(NEXTERNAL-1)
C     integer    n_max_cg
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'nsquaredSO.inc'
      character*7 pdlabel,epa_label
      integer lhaid
      common/to_pdf/lhaid,pdlabel,epa_label
C     
C     LOCAL
C     
      INTEGER I,J,K, iolp
C     four momenta. Energy is the zeroth component.
      REAL*8 P(0:3,NEXTERNAL-1)
      INTEGER MATELEM_ARRAY_DIM
      REAL*8 , ALLOCATABLE :: MATELEM(:,:)
      REAL*8 SQRTS,AO2PI,TOTMASS
C     sqrt(s)= center of mass energy 
      REAL*8 PIN(0:3), POUT(0:3)
      CHARACTER*140 BUFF(NEXTERNAL-1)
      CHARACTER*140 BUFF2
      CHARACTER*80 BUFF80
      INTEGER RETURNCODE, UNITS, TENS, HUNDREDS
      INTEGER NSQUAREDSO_LOOP
      REAL*8 , ALLOCATABLE :: PREC_FOUND(:)

C     
C     GLOBAL VARIABLES
C     
C     This is from ML code for the list of split orders selected by
C     the process definition
C     
      INTEGER NLOOPCHOSEN
      CHARACTER*20 CHOSEN_LOOP_SO_INDICES(NSQUAREDSO)
      LOGICAL CHOSEN_LOOP_SO_CONFIGS(NSQUAREDSO)
      COMMON/ML5_0_CHOSEN_LOOP_SQSO/CHOSEN_LOOP_SO_CONFIGS
      INTEGER NBORNCHOSEN
      CHARACTER*20 CHOSEN_BORN_SO_INDICES(NSQSO_BORN)
      LOGICAL CHOSEN_BORN_SO_CONFIGS(NSQSO_BORN)
      COMMON/ML5_0_CHOSEN_BORN_SQSO/CHOSEN_BORN_SO_CONFIGS


      double precision tolerance
      parameter (tolerance=1d-4)

      integer nwgt,max_rwgt
      parameter(max_rwgt=200)
      double precision rwgtxsec(max_rwgt)
      common /rwgt_wgts/rwgtxsec,nwgt


      integer ievent,ii,ifile,maxparticles,iord,cpower(3),iiwgt,iii
      data cpower/0,1,2/
      parameter (maxparticles=10)
      double precision virtual_Xsec(0:3,max_rwgt),virtual_unc(0:3)
     $     ,rwgt_muR_dep_fac,scale,wgt_to_add,alphas
      logical done
      INTEGER NUP,IDPRUP,IDUP(maxparticles),ISTUP(maxparticles),MOTHUP(2
     $     ,maxparticles),ICOLUP(2,maxparticles)
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP, PUP(5,maxparticles)
     $     ,VTIMUP(maxparticles),SPINUP(maxparticles)

      include 'c_weight.inc'

      integer ntot,nsun,nsps,nups,neps,n100,nddp,nqdp,nini,n10,n1(0:9)

C     
C     SAVED VARIABLES
C     
      LOGICAL INIT
      DATA INIT/.TRUE./
      COMMON/INITCHECKSA/INIT
C     
C     EXTERNAL
C     
      REAL*8 DOT
      EXTERNAL DOT

! mio 
! permutation crap
c      integer mom_permutation(nexternal)
c      data mom_permutation / 1, 2, 5, 3, 4, 6/
! end permutation crap`
! mio  
      double precision new_wgt, sum_wgt, sum_wgt_new, unc_sum_new
      integer ifilerwgt

C     
C     BEGIN CODE
C     
C     

      IF (INIT) THEN
        INIT=.FALSE.
        CALL ML5_0_GET_ANSWER_DIMENSION(MATELEM_ARRAY_DIM)
        ALLOCATE(MATELEM(0:3,0:MATELEM_ARRAY_DIM))
        CALL ML5_0_GET_NSQSO_LOOP(NSQUAREDSO_LOOP)
        ALLOCATE(PREC_FOUND(0:NSQUAREDSO_LOOP))

C       INITIALIZATION CALLS
C       
C       Call to initialize the values of the couplings, masses and
C        widths 
C       used in the evaluation of the matrix element. The primary
C        parameters of the
C       models are read from Cards/param_card.dat. The secondary
C        parameters are calculated
C       in Source/MODEL/couplings.f. The values are stored in common
C        blocks that are listed
C       in coupl.inc .
C       first call to setup the paramaters
        CALL SETPARA('../Cards/param_card.dat')
        pdlabel='lhapdf'
*       pdlabel='nn23nlo'
        lhaid=90500
        call pdfwrap
C       set up masses
        INCLUDE 'pmass.inc'

      ENDIF

      sum_wgt=0d0
      sum_wgt_new=0d0
      unc_sum_new=0d0


C     Start by initializing what is the squared split orders indices
C      chosen
      NLOOPCHOSEN=0
      DO I=1,NSQUAREDSO
        IF (CHOSEN_LOOP_SO_CONFIGS(I)) THEN
          NLOOPCHOSEN=NLOOPCHOSEN+1
          WRITE(CHOSEN_LOOP_SO_INDICES(NLOOPCHOSEN),'(I3,A2)') I,'L)'
        ENDIF
      ENDDO
      NBORNCHOSEN=0
      DO I=1,NSQSO_BORN
        IF (CHOSEN_BORN_SO_CONFIGS(I)) THEN
          NBORNCHOSEN=NBORNCHOSEN+1
          WRITE(CHOSEN_BORN_SO_INDICES(NBORNCHOSEN),'(I3,A2)') I,'B)'
        ENDIF
      ENDDO

      MadLoopResourcePath = 'MadLoop5_resources'
      CALL SETMADLOOPPATH(MadLoopResourcePath)

      CALL PRINTOUT()
      
      ifile=12
      open (unit=ifile,file='events.lhe',status='old')

      ifilerwgt=14
      open (unit=ifilerwgt,file='events_rwgt.lhe',status='unknown')

      do
         read(ifile,'(a)') buff80
         write(ifilerwgt,'(a)') buff80
         if (index(buff80,"</init>").ne.0) exit
      enddo

      ievent=0
      do iord=0,3
         do iiwgt=1,max_rwgt
            virtual_Xsec(iord,iiwgt)=0d0
         enddo
         virtual_unc(iord)=0d0
      enddo
      do
         call read_lhef_event(ifile,NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP
     $        ,AQCDUP,IDUP,ISTUP,MOTHUP,ICOLUP,PUP,VTIMUP,SPINUP,buff2
     $        ,done)
         if (done) exit

         ievent=ievent+1
         if (mod(ievent,100).eq.0) then
            write (*,*) 'at event',ievent,':',sum_wgt_new
     $           /dble(ievent),'+/-',sqrt(unc_sum_new)/dble(ievent)
            call flush
         endif
         call fill_wgt_info_from_rwgt_lines

         do i=1,icontr
            g=g_strong(i)
            MU_R=sqrt(scales2(1,i)) ! ellis-sexton scale
            call update_as_param()
            
            do j=1,nexternal-1
               do ii=0,3 ! mio
c! permutation crap
c                  p(ii,mom_permutation(j))=momenta_m(ii,j,1,i)
c! end permutation crap
                   p(ii,j)=momenta_m(ii,j,1,i)
               enddo
            enddo

      include 'check_olp.inc'

c Update the statistics using the MadLoop return code (returncode)
            ntot = ntot+1       ! total number of PS
            if (returncode/100.eq.1) then
               nsun = nsun+1    ! stability unknown
            elseif (returncode/100.eq.2) then
               nsps = nsps+1    ! stable PS point
            elseif (returncode/100.eq.3) then
               nups = nups+1    ! unstable PS point, but rescued
            elseif (returncode/100.eq.4) then
               neps = neps+1    ! exceptional PS point: unstable, and not possible to rescue
            else
               n100=n100+1      ! no known returncode (100)
            endif
            if (mod(returncode,100)/10.eq.1 .or. mod(returncode,100)
     $           /10.eq.3) then
               nddp = nddp+1    ! only double precision was used
               if (mod(returncode,100)/10.eq.3) nini=nini+1 ! MadLoop initialization phase
            elseif (mod(returncode,100)/10.eq.2 .or. mod(returncode,100)/10.eq.4)
     $              then
               nqdp = nqdp+1    ! quadruple precision was used
               if (mod(returncode,100)/10.eq.4) nini=nini+1 ! MadLoop initialization phase
            else
               n10=n10+1        ! no known returncode (10)
            endif
            n1(mod(returncode,10))=n1(mod(returncode,10))+1 ! unit ret code distribution

            if ((matelem(0,0)-wgt_ME_tree(1,i))/(matelem(0,0)
     $           +wgt_ME_tree(1,i)).gt.0.001) then
               write (*,*) 'Large difference in Born:',matelem(0,0)
     $              ,wgt_ME_tree(1,i)
            endif

C take the original event weight, multiply it by loop ME and
C divide by the born
            new_wgt=xwgtup*matelem(1,0)/matelem(0,0)
            sum_wgt=sum_wgt+xwgtup
            sum_wgt_new=sum_wgt_new+new_wgt
            unc_sum_new=unc_sum_new+new_wgt**2

            call write_lhef_event(ifilerwgt,
     #             NUP,IDPRUP,new_wgt,SCALUP,AQEDUP,AQCDUP,
     #             IDUP,ISTUP,MOTHUP,ICOLUP,PUP,VTIMUP,SPINUP,buff)

            matelem(1,0)=0d0
c            do iiwgt=1,nwgt
c               do iord=1,0,-1
c                  wgt_to_add=rwgtxsec(iiwgt)/matelem(0,0)
c                  scale=sqrt(scales2(2,i))
c                  if (iiwgt.eq.2 .or. iiwgt.eq.5 .or. iiwgt.eq.8) then
c                     ! twice renormalisation scale
c                     wgt_to_add=wgt_to_add/alphas(scale)
c                     scale=sqrt(scales2(2,i))*2d0
c                     wgt_to_add=wgt_to_add*alphas(scale)
c                  elseif (iiwgt.eq.3 .or. iiwgt.eq.6 .or. iiwgt.eq.9) then
c                     ! half renormalisation scale
c                     wgt_to_add=wgt_to_add/alphas(scale)
c                     scale=sqrt(scales2(2,i))/2d0
c                     wgt_to_add=wgt_to_add*alphas(scale)
c                  endif
c                  if (iiwgt.eq.1 .and. iord.ne.0) then
c                     matelem(1,iord)=matelem(1,iord)
c     $                    *rwgt_muR_dep_fac(sqrt(scales2(2,i)),scale
c     $                    ,dble(cpower(iord)))
c                     matelem(1,0)=matelem(1,0)+matelem(1,iord)
c                  endif
c                  wgt_to_add=wgt_to_add*matelem(1,iord)
c
c                  if (returncode/100.eq.4) then
c                      ! Exceptional unstable PS point that has not been rescued.
c                      ! Simply use the accumulated average...
c                     if (iiwgt.eq.1 .and. iord.eq.0) write (*,*)
c     $                    'CHECK_OLP: USE THE ACCUMULATED VALUE'
c                     if (ievent.ne.0) then
c                        wgt_to_add=virtual_Xsec(iord,iiwgt)/dble(ievent)
c                     else
c                        wgt_to_add=0d0
c                     endif
c                  endif
c
c                  virtual_Xsec(iord,iiwgt)=virtual_Xsec(iord,iiwgt)
c     $                 +wgt_to_add
c                  if (iiwgt.eq.1) then
c                     virtual_unc(iord)=virtual_unc(iord)+(wgt_to_add)**2
c                  endif
c               enddo
c            enddo

         enddo
      enddo

      write (*,*) 'virtual cross section=',virtual_Xsec(0,1)
     $     /dble(ievent),'+/-',sqrt(virtual_unc(0))/dble(ievent)

      open (unit=13,file='events.virtual',status='unknown')
      do iiwgt=1,nwgt
         write (13,*) (virtual_xsec(iord,iiwgt)/dble(ievent),iord=0,3)
      enddo
      write (13,*) ''
      write(13,*) (sqrt(virtual_unc(iord))/dble(ievent),iord=0,3)

      DEALLOCATE(MATELEM)
      DEALLOCATE(PREC_FOUND)


      if (ntot.ne.0) then
         write(*,*) "Satistics from MadLoop:"
         write(*,*)
     &        "  Total points tried:                              ",ntot
         write(*,*)
     &        "  Stability unknown:                               ",nsun
         write(*,*)
     &        "  Stable PS point:                                 ",nsps
         write(*,*)
     &        "  Unstable PS point (and rescued):                 ",nups
         write(*,*)
     &        "  Exceptional PS point (unstable and not rescued): ",neps
         write(*,*)
     &        "  Double precision used:                           ",nddp
         write(*,*)
     &        "  Quadruple precision used:                        ",nqdp
         write(*,*)
     &        "  Initialization phase-space points:               ",nini
         write(*,*)
     &        "  Unknown return code (100):                       ",n100
         write(*,*)
     &        "  Unknown return code (10):                        ",n10
         write(*,*)
     &        "  Unit return code distribution (1):               "
         do j=0,9
           if (n1(j).ne.0) then
              write(*,*) "#Unit ",j," = ",n1(j)
           endif
         enddo
      endif

      write(*,*) 'SUM OF ORIGINAL WEIGHTS:', sum_wgt
      write(*,*) 'SUM OF NEW WEIGHTS:', sum_wgt_new

      close(ifile)
     
      write(ifilerwgt, '(a)') "</LesHouchesEvents>"
      close(ifilerwgt)


      END

















c     This function defines the reweighting of the cross section to
c     include a muR-dependent pre-factor. Multiply by the muR-dependent
c     factor and devide by the muR-independent one.
c Note: This is implemented below for the Bottom Yukawa in the SM.
c       Change it to the factor you need to reweight.
      Double precision function rwgt_muR_dep_fac(scale,central,cpower_in)
      implicit none
      double precision cpower_in
      double precision scale,vev,mbmb,apimuR,apimZ,apimb,mbmuR,alphas,pi
      parameter (pi=3.14159265358979323846d0)
      include "nexternal.inc"
c      include "genps.inc"
c      include "reweight.inc"
      include "coupl.inc"
      include "../Source/MODEL/input.inc"
      integer i
      double precision central,tootiny,apicentral,mbcentral
      parameter (tootiny=1d-9)
      rwgt_muR_dep_fac = 1d0
c     This is relevant for a muR-dependent bottom-mass in Yukawa.
      IF(cpower_in .ne. 0d0) THEN
      vev    = 246.21845813429518469305d0 !vev in aMC@NLO from y_b->m_b
      mbmb = MDL_YB*vev/dsqrt(2d0)
com-- mbmb input for fixed Yukawa bmass in param_card.dat is used here
com-- as start value of running and to remove it from the cross section
c new settings NLO
      apimZ  = alphas(MDL_MZ)/pi

      if(dabs(scale/central-1d0).lt.tootiny) then
c if scale muR is the same as the central scale of muR, get 
c "input value" mb(muR) with highest possible accuracy
         CALL runalpha(apimZ,MDL_MZ,central,4d0,4,0,apimuR)
         CALL runalpha(apimZ,MDL_MZ,mbmb,4d0,4,0,apimb)
         CALL runmass(mbmb,apimb,apimuR,4d0,4,mbmuR)
      else
c if scale and central are different (muR variations) do two steps:
c step 1: get "input value" mb(central scale) from most accurate running
         CALL runalpha(apimZ,MDL_MZ,central,4d0,4,0,apicentral)
         CALL runalpha(apimZ,MDL_MZ,mbmb,4d0,4,0,apimb)
         CALL runmass(mbmb,apimb,apicentral,4d0,4,mbcentral)
c step 2: get variation around central value, ie mb(muR), with loop 
c         order consistent with computation LO: 1-loop, NLO: 2-loop
         CALL runalpha(apicentral,central,scale,4d0,2,0,apimuR)
         CALL runmass(mbcentral,apicentral,apimuR,4d0,2,mbmuR)
      endif
      rwgt_muR_dep_fac = (mbmuR/mbmb)**cpower_in
      ELSE
         return
      ENDIF
      END

C-{{{ routines for running of alphas:

C-{{{ subroutine odeint:

      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
c..(C) Copr. 1986-92 Numerical Recipes Software 5,".
c..   transscribed to real*8 by R. Harlander, Feb.2002
      implicit real*8 (a-z)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(dabs(x-xsav).gt.dabs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(dabs(hnext).lt.hmin) write(6,*)
     *       'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      write(6,*) 'too many steps in odeint'
      stop
      return
      END

C-}}}
C-{{{ subroutine rkck:

      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
c..(C) Copr. 1986-92 Numerical Recipes Software 5,".
c..   transscribed to real*8 by R. Harlander, Feb.2002
      implicit real*8 (a-z)
      INTEGER n,NMAX
      REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs
      INTEGER i
      REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END

C-}}}
C-{{{ subroutine rkqs:

      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
c..(C) Copr. 1986-92 Numerical Recipes Software 5,".
c..   transscribed to real*8 by R. Harlander, Feb.2002
      implicit real*8 (a-z)
      INTEGER n,NMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER i
      REAL*8 errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     *PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89d-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.d0
      do 11 i=1,n
        errmax=max(errmax,dabs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(dabs(htemp),0.1*dabs(h)),h)
        xnew=x+h
        if(xnew.eq.x) write(6,*) 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END

C-}}}
C-{{{ subroutine runalpha:

      subroutine runalpha(api0,mu0,mu,nf,nloop,verb,apiout)
C..
c..   NEEDS:  rkck.f rkqs.f odeint.f  (from Numerical Recipes)
c..   
c..   Note:  api = {\alpha_s \over \pi}
C..
c..   purpose : computes the value of api(mu) from api(mu0)
c..   method  : solving RG-equation by adaptive Runge-Kutta method
c..   uses    : odeint.for  from Numerical Recipes
C..
c..   api0  :  api(mu0)
c..   nf    :  number of flavors
c..   nloop :  number of loops
c..   verb  :  0=quiet,  1=verbose
c..   apiout:  api(mu)    
C..
      implicit real*8 (a-h,o-z)
      INTEGER KMAXX,NMAX,NVAR
      PARAMETER (KMAXX=200,NMAX=50,NVAR=1)
      INTEGER kmax,kount,nbad,nok,nrhs,verb
      REAL*8 dxsav,eps,h1,hmin,x,y,apif(NVAR),api0,apiout,pi
      real*8 mu,mu0,l0,lf,nf
c..   /path/  is for odeint.for:
      COMMON /path/ kmax,kount,dxsav,x(KMAXX),y(NMAX,KMAXX)
      common /bfunc/ beta0,beta1,beta2,beta3
      COMMON /cbnrhs/nrhs 
      data pi/3.14159265358979323846264338328d0/
      EXTERNAL rhs,rkqs

      if (nloop.eq.0) then
         apiout = api0
         return
      endif

      nrhs=0

c..   integration bounds (note that log(mu^2) is the integration variable)
      l0 = 0.d0
      lf = 2.*dlog(mu/mu0)
      apif(1)=api0

c..   see documentation for odeint.for:
      eps=1.0d-8
      h1=dabs(lf-l0)/10.d0
      hmin=0.0d0
      kmax=100
      dxsav=dabs(lf-l0)/20.d0

c..   initialize beta-function (common block /bfunc/):
      call inibeta(nf,nloop)

c..   check if input values are reasonable
      dlam = mu0*dexp(-1.d0/(2.d0*beta0*api0))
      if (mu.le.dlam) then
         write(6,2001) dlam,mu,mu0,api0*pi
      endif

c..   integrate RG-equation:

      call odeint(apif,NVAR,l0,lf,eps,h1,hmin,nok,nbad,rhs,rkqs)

      if (verb.eq.1) then
         write(6,'(/1x,a,t30,i3)') 'Successful steps:',nok
         write(6,'(1x,a,t30,i3)') 'Bad steps:',nbad
         write(6,'(1x,a,t30,i3)') 'Function evaluations:',nrhs
         write(6,'(1x,a,t30,i3)') 'Stored intermediate values:',kount
      endif

c..   api(mu):
      apiout = apif(1)

 2001 format(' -> <subroutine runalpha>',/,
     &     ' - WARNING: mu-value too low.',/,
     &     ' -     Should be significantly larger than  ',1f8.3,'.',/,
     &     ' -             mu = ',1f8.3,' GeV',/,
     &     ' -            mu0 = ',1f8.3,' GeV',/,
     &     ' -        api0*pi = ',1f8.3,/,
     &     ' -     Integration might break down.',/,
     &     '<- <subroutine runalpha>'
     &     )

      END

C-}}}
C-{{{ subroutine rhs:

      subroutine rhs(logmumu0,api,ainteg)
C..
c..   RG-equation:   (d api)/(d log(mu^2)) = api*beta(api)
C..
      implicit real*8 (a-h,o-z)
      integer nrhs
      real*8 api(*),ainteg(*),logmumu0
      common /bfunc/ beta0,beta1,beta2,beta3
      COMMON nrhs
      nrhs=nrhs+1
      ainteg(1) = api(1)*(- beta0*api(1) - beta1*api(1)**2 - beta2
     &     *api(1)**3 - beta3*api(1)**4)
      end

C-}}}
C-{{{ subroutine inibeta:

      subroutine inibeta(nf,nloopin)
C..
c..   initialize beta function
C..
      implicit real*8 (a-h,o-z)
      real*8 nf
      data z3/1.2020569031595942853997/
      common /bfunc/ beta0,beta1,beta2,beta3

      beta0 = (33 - 2*nf)/12.d0
      beta1 = (102 - (38*nf)/3.d0)/16.d0
      beta2 = (2857/2.d0 - (5033*nf)/18.d0 + (325*nf**2)/54.d0)/64.d0
      beta3 = (149753/6.d0 + (1093*nf**3)/729.d0 + 3564*z3 + nf**2
     &     *(50065/162.d0 + (6472*z3)/81.d0) - nf*(1078361/162.d0 +
     &     (6508*z3)/27.d0))/256.d0

      
      nloop=nloopin

      if (nloop.gt.4) then
         write(6,*) '-> <subroutine inibeta>:'
         write(6,*)
     &        ' - 5-loop beta function unknown. Using 4-loop instead.'
         write(6,*) '<- <subroutine inibeta>'
         nloop=4
      endif
      if (nloop.lt.4) then
         beta3 = 0d0
         if (nloop.lt.3) then
            beta2 = 0d0
            if (nloop.lt.2) then
               beta1 = 0d0
               if (nloop.lt.1) then
                  beta0=0d0
               endif
            endif
         endif
      endif
      end

C-}}}

C-}}}
C-{{{ subroutine runmass:

      subroutine runmass(mass0,api0,apif,nf,nloop,massout)
c..
c..   evaluates the running of the MS-bar quark mass
c..   by expanding the equation
c..   
c..   m(mu) = m(mu0) * exp( \int_a0^af dx gammam(x)/x/beta(x) )
c..   
c..   in terms of alpha_s. The results agree with RunDec.m.
c..   
c..   
c..   Input:
c..   ------
c..   mass0  :  m(mu0)
c..   api0   :  alpha_s(mu0)/pi
c..   apif   :  alpha_s(muf)/pi
c..   nf     :  number of flavors
c..   nloop  :  order of calculation (nloop=1..4)
c..
c..   Output:
c..   -------
c..   massout:  m(muf)
c..   
      implicit real*8 (a-h,o-z)
      real*8 mass0,massout,massfun
      real*8 nf
      external massfun
      parameter(accmass=1.d-6)
      common /bfunc/ beta0,beta1,beta2,beta3
      common /gfunc/ gamma0,gamma1,gamma2,gamma3

      if (nloop.eq.0) then
         massout = mass0
         return
      endif

      call inigamma(nf,nloop)
      call inibeta(nf,nloop)

      bb1 = beta1/beta0
      bb2 = beta2/beta0
      bb3 = beta3/beta0

      cc0 = gamma0/beta0
      cc1 = gamma1/beta0
      cc2 = gamma2/beta0
      cc3 = gamma3/beta0

      cfunc1 = 1.d0
      cfunc2 = cc1 - bb1*cc0
      cfunc3 = 1/2.d0*((cc1-bb1*cc0)**2 + cc2 - bb1*cc1 + bb1**2*cc0 -
     &     bb2*cc0)
      cfunc4 = (1/6*(cc1 - bb1*cc0)**3 + 1/2*(cc1 - bb1*cc0)*(cc2 - bb1
     &     *cc1 + bb1**2*cc0 - bb2*cc0) + 1/3*(cc3 - bb1*cc2 + bb1**2
     &     *cc1 - bb2*cc1 - bb1**3*cc0 + 2*bb1*bb2*cc0 - bb3*cc0))

      if (nloop.lt.4) then
         cfunc4 = 0.d0
         if (nloop.lt.3) then
            cfunc3 = 0.d0
            if (nloop.lt.2) then
               cfunc2 = 0.d0
               if (nloop.lt.1) then
                  cfunc1 = 0.d0
               endif
            endif
         endif
      endif

      cfuncmu0 = cfunc1 + cfunc2*api0 + cfunc3*api0**2 + cfunc4*api0**3
      cfuncmuf = cfunc1 + cfunc2*apif + cfunc3*apif**2 + cfunc4*apif**3


      massout = mass0*(apif/api0)**cc0*cfuncmuf/cfuncmu0
      
      return
      end

C-}}}
C-{{{ subroutine inigamma:

      subroutine inigamma(nfin,nloopin)
C
C     initialize beta function
C
      implicit real*8 (a-h,o-z)
      real*8 nf,nfin
      data z3/1.2020569031595942853997/,
     &     z5/1.0369277551433699263/,
     &     pi/3.1415926535897932381/
      common /gfunc/ gamma0,gamma1,gamma2,gamma3

      nf = nfin

      gamma0 = 1.d0
      gamma1 = (67.33333333333333d0 - (20*nf)/9.d0)/16.d0
      gamma2 = (1249.d0 - (140*nf**2)/81.d0 + 2*nf*(-20.59259259259259d0
     &     - 48*z3) +(8*nf*(-46 + 48*z3))/9.d0)/64.d0
      gamma3 = (28413.91975308642d0 + (135680*z3)/27.d0 + nf**3*(-1
     &     .3662551440329218d0 + (64*z3)/27.d0) + nf**2*(21
     &     .57201646090535d0 - (16*Pi**4)/27.d0 + (800*z3)/9.d0) - 8800
     &     *z5 + nf*(-3397.1481481481483d0 + (88*Pi**4)/9.d0 - (34192
     &     *z3)/9.d0 + (18400*z5)/9.d0))/256.d0

      nloop=nloopin

      if (nloop.gt.4) then
         write(6,*) '-> <subroutine inigamma>:'
         write(6,*)
     &        ' - 5-loop gamma function unknown. Using 4-loop instead.'
         write(6,*) '<- <subroutine inigamma>'
         nloop=4
      endif
      if (nloop.lt.4) then
         gamma3 = 0d0
         if (nloop.lt.3) then
            gamma2 = 0d0
            if (nloop.lt.2) then
               gamma1 = 0d0
               if (nloop.lt.1) then
                  gamma0 = 0d0
               endif
            endif
         endif
      endif
      end


C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION ALPHAS(Q)
C     wrapper to the lhapdf alphaS
C-----------------------------------------------------------------------------
      IMPLICIT NONE
c
      REAL*8 Q,alphasPDF
      external alphasPDF

c timing statistics

      ALPHAS=alphasPDF(Q)


      RETURN
      END


      subroutine write_lhef_event(ifile,
     # NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     # IDUP,ISTUP,MOTHUP,ICOLUP,PUP,VTIMUP,SPINUP,buff)
      implicit none
      INTEGER NUP,IDPRUP,IDUP(*),ISTUP(*),MOTHUP(2,*),ICOLUP(2,*)
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,
     # PUP(5,*),VTIMUP(*),SPINUP(*)
      character*140 buff
      integer ifile,i,kk
      character*9 ch1
      integer isorh_lhe,ifks_lhe,jfks_lhe,fksfather_lhe,ipartner_lhe
      double precision scale1_lhe,scale2_lhe
      integer ii,j,nps,nng,iFKS,idwgt
      double precision wgtcentral,wgtmumin,wgtmumax,wgtpdfmin,wgtpdfmax
      integer event_id
      common /c_event_id/ event_id
      integer i_process
      common/c_i_process/i_process
      integer nattr,npNLO,npLO
      common/event_attributes/nattr,npNLO,npLO
c      include 'reweight_all.inc'
c      include './run.inc'
c      include 'unlops.inc'
c     if event_id is zero or positive (that means that there was a call
c     to write_lhef_header_banner) update it and write it
c RF: don't use the event_id:
      event_id = -99
c
c      if (event_id.ge.0) then
c         event_id=event_id+1
c         if (event_id.le.9) then
c            write(ifile,'(a,i1,a)') "  <event id='",event_id,"'>"
c         elseif(event_id.le.99) then
c            write(ifile,'(a,i2,a)') "  <event id='",event_id,"'>"
c         elseif(event_id.le.999) then
c            write(ifile,'(a,i3,a)') "  <event id='",event_id,"'>"
c         elseif(event_id.le.9999) then
c            write(ifile,'(a,i4,a)') "  <event id='",event_id,"'>"
c         elseif(event_id.le.99999) then
c            write(ifile,'(a,i5,a)') "  <event id='",event_id,"'>"
c         elseif(event_id.le.999999) then
c            write(ifile,'(a,i6,a)') "  <event id='",event_id,"'>"
c         elseif(event_id.le.9999999) then
c            write(ifile,'(a,i7,a)') "  <event id='",event_id,"'>"
c         elseif(event_id.le.99999999) then
c            write(ifile,'(a,i8,a)') "  <event id='",event_id,"'>"
c         elseif(event_id.le.999999999) then
c            write(ifile,'(a,i9,a)') "  <event id='",event_id,"'>"
c         else
c            write (ifile,*) "ERROR: EVENT ID TOO LARGE",event_id
c            write (*,*) "ERROR: EVENT ID TOO LARGE",event_id
c            stop
c         endif
c      elseif(nattr.eq.2) then
c         if ( (npLO.ge.10.or.npLO.lt.0) .and.
c     &        (npNLO.ge.10.or.npNLO.lt.0)) then
c            write(ifile,'(a,i2,a,i2,a)') "  <event npLO=' ",npLO
c     $           ," ' npNLO=' ",npNLO," '>"
c         elseif( (npLO.lt.10.or.npLO.ge.0) .and.
c     &        (npNLO.ge.10.or.npNLO.lt.0)) then
c            write(ifile,'(a,i1,a,i2,a)') "  <event npLO=' ",npLO
c     $           ," ' npNLO=' ",npNLO," '>"
c         elseif( (npLO.ge.10.or.npLO.lt.0) .and.
c     &        (npNLO.lt.10.or.npNLO.ge.0)) then
c            write(ifile,'(a,i2,a,i1,a)') "  <event npLO=' ",npLO
c     $           ," ' npNLO=' ",npNLO," '>"
c         elseif( (npLO.lt.10.or.npLO.ge.0) .and.
c     &        (npNLO.lt.10.or.npNLO.ge.0)) then
c            write(ifile,'(a,i1,a,i1,a)') "  <event npLO=' ",npLO
c     $           ," ' npNLO=' ",npNLO," '>"
c         endif
c      else
         write(ifile,'(a)') '  <event>'
c      endif
      write(ifile,503)NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
      do i=1,nup
        write(ifile,504)IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I),
     #                  ICOLUP(1,I),ICOLUP(2,I),
     #                  PUP(1,I),PUP(2,I),PUP(3,I),PUP(4,I),PUP(5,I),
     #                  VTIMUP(I),SPINUP(I)
      enddo
c      if(buff(1:1).eq.'#' .and. (do_rwgt_scale .or. do_rwgt_pdf .or.
c     &     jwgtinfo.lt.0)) then
c         write(ifile,'(a)') buff(1:len_trim(buff))
c         read(buff,*)ch1,iSorH_lhe,ifks_lhe,jfks_lhe,
c     #                    fksfather_lhe,ipartner_lhe,
c     #                    scale1_lhe,scale2_lhe,
c     #                    jwgtinfo,mexternal,iwgtnumpartn,
c     #         wgtcentral,wgtmumin,wgtmumax,wgtpdfmin,wgtpdfmax
c         if(jwgtinfo.eq.-5.or.jwgtinfo.eq.-9) then
c            write(ifile,'(a)')'  <mgrwgt>'
c            write (ifile,'(1x,d16.10,3(1x,i4))') wgtref,n_ctr_found
c     &           ,n_mom_conf, nint(wgtcpower)
c            do i=1,n_mom_conf
c               do j=1,mexternal
c                  write (ifile,'(4(1x,d21.15))')
c     &                 (momenta_str(ii,j,i),ii=0,3)
c               enddo
c            enddo
c            do i=1,n_ctr_found
c               write (ifile,'(a)') trim(adjustl(n_ctr_str(i)))
c            enddo
c            write(ifile,'(a)')'  </mgrwgt>'
c         endif
c         if(jwgtinfo.eq.15) then
c            write(ifile,'(a)')'  <unlops>'
c            write(ifile,*)NUP_H
c            do i=1,NUP_H
c               write(ifile,504)IDUP_H(I),ISTUP_H(I),MOTHUP_H(1,I)
c     $              ,MOTHUP_H(2,I),ICOLUP_H(1,I),ICOLUP_H(2,I),PUP_H(1
c     $              ,I),PUP_H(2,I),PUP_H(3,I),PUP_H(4,I),PUP_H(5,I),
c     $              VTIMUP_H(I),SPINUP_H(I)
c            enddo
c            write(ifile,'(a)')'  </unlops>'
c         endif
c         if(abs(jwgtinfo).eq.9)then
c            if (do_rwgt_scale .or. do_rwgt_pdf) then
c               write(ifile,'(a)') '  <rwgt>'
c               idwgt=1000
c               if (do_rwgt_scale) then
c                  do kk=1,dyn_scale(0)
c                     if (lscalevar(kk)) then
c                        do i=1,nint(scalevarF(0))
c                           do j=1,nint(scalevarR(0))
c                              idwgt=idwgt+1
c                              write(ifile,601) "   <wgt id='",idwgt,"'>"
c     $                             ,wgtxsecmu(j,i,kk)," </wgt>"
c                           enddo
c                        enddo
c                     else
c                        idwgt=idwgt+1
c                        write(ifile,601) "   <wgt id='",idwgt,"'>"
c     $                       ,wgtxsecmu(1,1,kk)," </wgt>"
c                     endif
c                  enddo
c               endif
c               if (do_rwgt_pdf) then
c                  do j=1,lhaPDFid(0)
c                     if (lpdfvar(j)) then
c                        do i=0,nmemPDF(j)
c                           idwgt=idwgt+1
c                           write(ifile,601) "   <wgt id='",idwgt,"'>"
c     $                          ,wgtxsecPDF(i,j)," </wgt>"
c                        enddo
c                     else
c                        idwgt=idwgt+1
c                        write(ifile,601) "   <wgt id='",idwgt,"'>"
c     $                       ,wgtxsecPDF(0,j)," </wgt>"
c                     endif
c                  enddo
c               endif
c               write(ifile,'(a)') '  </rwgt>'
c            endif
c         endif
c      endif
      write(ifile,'(a)') '  </event>'
 401  format(2(1x,e14.8))
 402  format(8(1x,e14.8))
 403  format(6(1x,e14.8))
 404  format(3(1x,e14.8))
 405  format(4(1x,e14.8))
 406  format(2(1x,e14.8),2(1x,i3))
 441  format(4(1x,e16.10))
 442  format(1x,e16.10,2(1x,e14.8))
 503  format(1x,i2,1x,i6,4(1x,e14.8))
 504  format(1x,i8,1x,i2,4(1x,i4),5(1x,e14.8),2(1x,e10.4))
 601  format(a12,i4,a2,1x,e11.5,a7)
c
      return
      end






      subroutine read_lhef_event(ifile,
     # NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     # IDUP,ISTUP,MOTHUP,ICOLUP,PUP,VTIMUP,SPINUP,buff,done)
      implicit none
      logical done
      INTEGER NUP,IDPRUP,IDUP(*),ISTUP(*),MOTHUP(2,*),ICOLUP(2,*)
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,
     # PUP(5,*),VTIMUP(*),SPINUP(*)
      integer ifile,i,kk
      character*140 buff
      character*80 string
      character*9 ch1
      integer isorh_lhe,ifks_lhe,jfks_lhe,fksfather_lhe,ipartner_lhe
     $     ,mexternal,iwgtnumpartn,jwgtinfo
      double precision scale1_lhe,scale2_lhe
      integer ii,j,idwgt
      double precision wgtcentral,wgtmumin,wgtmumax,wgtpdfmin,wgtpdfmax
c     
      double precision wgtref
      integer wgtcpower
      common/event_info/wgtref,wgtcpower

      integer n_ctr_found,max_n_ctr,n_mom_conf,maxparticles
      parameter (max_n_ctr=256,maxparticles=10)
      character*1024 n_ctr_str(max_n_ctr)
      double precision momenta_str(0:3,maxparticles,max_n_ctr)
      common /c_rwgt_lines/n_ctr_str,momenta_str,n_ctr_found,n_mom_conf
c
      integer nwgt,max_rwgt
      parameter(max_rwgt=200)
      double precision rwgtxsec(max_rwgt)
      common /rwgt_wgts/rwgtxsec,nwgt
c
      logical done2
c
      done=.false.
      do 
         read(ifile,'(a)',err=99,end=99)string
         if (index(string,'<event>').ne.0) exit
      enddo

      read(ifile,*)NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
      do i=1,nup
        read(ifile,*)IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I),
     #                  ICOLUP(1,I),ICOLUP(2,I),
     #                  PUP(1,I),PUP(2,I),PUP(3,I),PUP(4,I),PUP(5,I),
     #                  VTIMUP(I),SPINUP(I)
      enddo
      read(ifile,'(a)')buff
      read(buff,*)ch1,iSorH_lhe,ifks_lhe,jfks_lhe,
     #                 fksfather_lhe,ipartner_lhe,
     #                 scale1_lhe,scale2_lhe,
     #                 jwgtinfo,mexternal,iwgtnumpartn,
     #      wgtcentral,wgtmumin,wgtmumax,wgtpdfmin,wgtpdfmax
      if(jwgtinfo.eq.-5 .or. jwgtinfo.eq.-9) then
         read(ifile,'(a)')string
         read(ifile,*) wgtref,n_ctr_found,n_mom_conf,wgtcpower
         do i=1,n_mom_conf
            do j=1,mexternal
               read (ifile,*) (momenta_str(ii,j,i),ii=0,3)
            enddo
         enddo
         do i=1,n_ctr_found
            read (ifile,'(a)') n_ctr_str(i)
         enddo
         read(ifile,'(a)')string
      endif
      if(abs(jwgtinfo).eq.9)then
         read(ifile,'(a)')string
         nwgt=0
         do
            nwgt=nwgt+1
            call read_rwgt_line(ifile,idwgt,rwgtxsec(nwgt),done2)
            if (done2) then
               nwgt=nwgt-1
               exit
            endif
         enddo
      endif
      do 
         read(ifile,'(a)')string
         if (index(string,'</event').ne.0) exit
      enddo
c
      return
 99   done=.true.
      end



      subroutine fill_wgt_info_from_rwgt_lines
      implicit none
      include 'nexternal.inc'
      include 'c_weight.inc'
      integer i,idum,j,k,momenta_conf(2),ii
      integer n_ctr_found,max_n_ctr,n_mom_conf,maxparticles
      parameter (max_n_ctr=256,maxparticles=10)
      character*1024 n_ctr_str(max_n_ctr)
      double precision momenta_str(0:3,maxparticles,max_n_ctr)
      common /c_rwgt_lines/n_ctr_str,momenta_str,n_ctr_found,n_mom_conf
      icontr=n_ctr_found
      iwgt=1
      do i=1,icontr
         read(n_ctr_str(i),*)(wgt(j,i),j=1,3),(wgt_ME_tree(j,i),j=1,2)
     &        ,idum,(pdg(j,i),j=1,nexternal),QCDpower(i),(bjx(j,i),j=1
     &        ,2),(scales2(j,i),j=1,3),g_strong(i),(momenta_conf(j),j=1
     &        ,2),itype(i),nFKS(i),idum,idum,idum,wgts(1,i)
         do ii=1,2
            do j=1,nexternal
               do k=0,3
                  if (momenta_conf(ii).gt.0) then
                     momenta_m(k,j,ii,i)=momenta_str(k,j
     $                                               ,momenta_conf(ii))
                  else
                     momenta_m(k,j,ii,i)=-99d0
                     exit
                  endif
               enddo
            enddo
         enddo
      enddo
      end
      


      subroutine read_rwgt_line(unit,id,wgt,done)
c read a line in the <rwgt> tag. The syntax should be
c  <wgt id='1001'> 0.1234567e+01 </wgt>
c The id should be exactly 4 digits long.
      implicit none
      logical done
      integer unit,id,wgt_start,id_start
      double precision wgt
      character*140 buff
      done=.false.
      read (unit,'(a)') buff
c Use char() to make sure that the non-standard characters are compiler
c independent (char(62)=">", char(61)="=", char(39)="'")
      wgt_start=index(buff,CHAR(39)//CHAR(62))+2
      id_start=index(buff,'id'//CHAR(61)//CHAR(39))+4
      read (buff(id_start:100),'(i4)',err=99,end=99) id
      read (buff(wgt_start:100),*,err=99,end=99) wgt
      return
 99   done=.true.
      end

