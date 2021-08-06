ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c main MD part of the code with various subtleties (a big bug can
c occur in every line, be carefull!)
c the variable it determines the mode:
c                          -1 : M equilibration
c                           0 : M sampling
c                          >0 : M+ fragmentation       (most important)
c                        9999 : M+ fragmentation trials(not used anymore)
c                          -2 : Ar bombarding
c
c mdok is set true if the result should be taken
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine md(it,isec,nuc,nmax,pos,iat,mass,imass,
     .              chrg,grad,velo,velof,
     .              list,tstep,k,totdump,nfragexit,fragm,fragf,
     .              fragat,dumpstep,etempin,prog,mdok,achrg,aspin,axyz,
     .              Tsoll,eimp,tadd,restart,Tav,Epav,Ekav,ttime,aTlast,
     .              fragstate,tfrag,method,ECP)
      implicit none

      integer nuc,iat(nuc),list(nuc),method
      integer nmax,it,k,totdump,prog,isec
      integer chrg,dumpstep,nfragexit,fragstate
      integer imass(nuc)
      real*8 pos (3,nuc)
      real*8 grad(3,nuc)
      real*8 velo(3,nuc)
      real*8 axyz(3,nuc)
      real*8 velof(nuc)
      real*8 aspin(nuc)
      real*8 achrg(nuc)
      real*8 mass (nuc)
      real*8 fragm(10)
      real*8 fragT(10)
      real*8 tstep,tadd,eimp,aTlast,tfrag
      real*8 T,Tav,Epav,Ekav,etempin,Tsoll,ttime
c 100 elements+100 for isotopes, 10 fragments max per fragmentation
      integer fragat(200,10)
      character*80 fragf(10)
      logical mdok,restart

      integer nstep,ndump,mdump,nfrag,i,j,avdump,kdump
      integer screendump,nadd,morestep,more,spin,fconst
      real*8 Epot,Eerror,Ekin,Edum,dum,etemp,Eav,fadd
      real*8 fstoau,angtoau,Ekinstart
      real*8 avspin(nuc),avchrg(nuc),avxyz(3,nuc)
      parameter (angtoau=1.0d0/0.529177260d0)
      parameter (fstoau =0.413413733365614D+02)
      character*2  asym
      character*20 fname
      logical err1,err2
      logical ECP
      include 'omp_lib.h'
      integer function omp_get_thread_num()
      integer thread
      thread = omp_get_thread_num()
      write(*,'(/13x''E N T E R I N G   M D   M O D U L E'',/)')

      mdok =.false.
      write(*,*)'mdok',mdok,thread
      nfrag=1
c status of fragments when run was ok (=0 is undefined)
c 1: normal
c 2: nothing happend for nfrag=2 for some time
      fragstate=0

c it=-1 : equilibration GS, no dump
c it= 0 : GS for generating starting points
c it> 1 : frag. MD

c MOLDEN file
      if(it.eq.0) then
                          fname='trjM'
      elseif(it.gt.0.and.it.lt.9999)then
                      write(fname,'(''trj.'',i4,''.'',i1)')it,isec
        if(it.lt.1000)write(fname,'(''trj.'',i3,''.'',i1)')it,isec
        if(it.lt.100) write(fname,'(''trj.'',i2,''.'',i1)')it,isec
        if(it.lt.10)  write(fname,'(''trj.'',i1,''.'',i1)')it,isec
      endif
      if(it.ge.0.and.it.lt.9999)open(unit=1,file=fname)

c this file is used to get the starting points = snapshots
      if(it.eq.0)then
         open(unit=2,file='qceims.gs')
         write(2,*)nmax
      endif
      write(*,*)'call ekinet',thread
c ini Ekin
      call ekinet(nuc,iat,velo,mass,Ekin,T)

      if(it.eq.9999)Tsoll=T
      Ekinstart=Ekin
c in frag runs the electronic temp is set
      if(it.gt.0.and.it.lt.9999)then
        write(*,*)'call setetemp2',thread
         call setetemp2(nfrag,eimp,etemp)
      else
         etemp=etempin
      endif
c CID mopac-pm6 behaves better without huge Electronic Temp.
      if(prog.eq.1.and.method.eq.3)then
         write(*,*) 'NOTE:In MOPAC CID runs the electronic T. is set
     .   to constant 300 K'
         etemp=300.0d0
      endif


c ini Epot
      spin=0
      write(*,*)'call egrad',thread
      call egrad(.true.,nuc,pos,iat,chrg,spin,prog,
     .           etemp,Epot,grad,achrg,aspin,method,
     .     ECP)
      if(Epot.eq.0) return

c do more more steps in a fragmentation run when nfragexit frag.
c are already there (i.e. after forming 2 frags, a third often occurs fast)
      more    =250
      avdump  =50
      screendump=100
      Tav     =0
      Epav    =0
      Ekav    =0
      Edum    =0
      Eerror  =0
      nstep   =0
       avchrg  =0
      avspin  =0
      aTlast  =0
      morestep=0
      fconst  =0
      tfrag   =0

      ndump=dumpstep
      kdump=avdump
c no additional calcs in ftemp (trial) runs
      if(it.eq.9999)morestep=more+1

c add the e-hit energy in the first MD steps linearly
      if(it.gt.0.and.it.lt.9999)then
         write(*,'(''Eimp (eV) = '',F6.1,5x,
     .             ''tauIC (fs) = '',F6.0,6x,
     .             ''nstep = '',i7,/)')
     .   eimp*27.212,tadd/fstoau,nmax
         nadd=(tadd+tstep)/tstep-1
         fadd=tstep/(tadd+tstep)
         write(*,'(''avcycle = '',i4,''    more = '',i4,/)')avdump,more
      else
         nadd =0
         velof=1.0d0
         screendump=500
      endif
      mdump=screendump


      write(*,'(''step   time [fs]'',4x,''Epot'',7x,''Ekin'',7x,
     .     ''Etot'',4x,''error'',2x,''#F   eTemp   frag. T'')')


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c MD loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do k=1,nmax
        write(*,*)'call do k',k,thread
      T=Ekin/(0.5*3*nuc*0.316681534524639E-05)
c in frag runs the electronic temp is set
      if(it.le.0.or.it.eq.9999)then
         etemp=etempin
c        in mopac cid the elc temp should be 300
         if(prog.eq.1.and.method.eq.3)then
            etemp=300.0d0
         endif
      endif

      Tav =Tav+T
      if(nfrag.eq.1)fragT(1)=Tav/k
      Epav=Epav+Epot
      Ekav=Ekav+Ekin
      if(nstep.gt.nadd)then
         Edum=Edum+Epot+Ekin
         Eav =Edum/float(nstep-nadd)
      else
         Eav=Epot+Ekin
      endif
c check av. Etot
      Eerror=Eav-Epot-Ekin
      if(it.eq.9999.or.it.lt.0)Eerror=0
c an error ocurred (normally failed to achieve SCF)
      err1=Epot.eq.0
      err2=abs(Eerror).gt.0.1.and.it.ne.9999
      if(err1.or.err2) then
c in case of errors take the traj anyway if fragments have been produced
         if( (nfrag.gt.1.and.nfrag.le.4). or. isec.gt.1 ) then
           write(*,8000)
     .     nstep,nstep*tstep/fstoau,Epot,Ekin,Epot+Ekin,
     .     Eerror,nfrag,etemp,fragT(1:nfrag)
           write(*,*)'EXIT DUE TO SCF CONVERGENCE PROBLEM'
           write(*,*)'OR LARGE MD INTEGRATION ERROR.'
           write(*,*)'because already nfrag > 1 result is taken'
           write(*,*)
     .     '(error often occurs for large inter-fragment distances)'
           write(*,*)
     .     ' and is therefore not very meaningfull)'
            mdok=.true.
         else
            mdok=.false.
            write(*,*) 'E,error ',Epot,Eerror
         endif

         goto 1000
      endif
      write(*,*)'end if',thread
c average over last avdump steps
      if(kdump.gt.avdump-1)then
         kdump =0
         avspin=0
         avchrg=0
         avxyz =0
         aTlast=0
      endif
      avchrg = avchrg + achrg
      avspin = avspin + aspin
      avxyz  = avxyz  + pos
      aTlast = aTlast + T

c print out every screendump steps
      if(mdump.gt.screendump-1)then
         write(*,8000)
     .   nstep,ttime,Epot,Ekin,Epot+Ekin,
     .   Eerror,nfrag,etemp,fragT(1:nfrag)
         mdump=0
      endif

c dump for MOLDEN
      if(ndump.gt.dumpstep-1.and.it.ge.0.and.it.lt.9999)then
         ndump=0
         write(1,*)nuc
         write(1,*)Epot
         do i=1,nuc
         write(1,'(1x,a2,1x,3F14.6)')
     .   asym(iat(i)),(pos(j,i)/angtoau,j=1,3)
         enddo
      endif

c no fragment run, dump to qceims.gs
      if(it.eq.0) then
         totdump=totdump+1
         do i=1,nuc
            write(2,'(6d16.8)')(pos(j,i),j=1,3),(velo(j,i),j=1,3)
         enddo
      endif

c do MD step
      call leapfrog(nuc,iat,grad,mass,tstep,pos,velo,Ekin,nstep)
c total time including secondaries
      ttime=ttime+tstep/fstoau
c calc E and forces
      call egrad(.false.,nuc,pos,iat,chrg,spin,prog,
     .           etemp,Epot,grad,achrg,aspin,method,ECP)

      ndump=ndump+1
      mdump=mdump+1
      kdump=kdump+1

c rescale to get Tsoll (NVT ensemble) for equilibration
c the GS sampling is done in NVE
      dum=100.*abs(Tav/k-Tsoll)/Tsoll
      if(dum.gt.5.0d0.and.it.lt.0.and.k.gt.50)then
         velo = velo/sqrt(Tav/k/Tsoll)
      endif

c add the IEE but only if not already fragmented
      if(it.gt.0.and.nstep.le.nadd.and.nfrag.eq.1)then
         call impactscsale(nuc,iat,velo,mass,velof,eimp,
     .                    fadd*nstep,Ekinstart)
      endif
c reduce eTemp when system is heated up (but only for nfrag=1)
      dum=eimp-eimp*float(nstep)/float(nadd)
      call setetemp2(nfrag,dum,etemp)
c Etemp for mopac is not performing well! Set it to fixed 300
      if(prog.eq.1.and.method.eq.3)then
         etemp=300.0d0
      endif
c increase temp til it breaks somewhere
      if(it.eq.9999)then
         Tsoll= Tsoll + 15.0
         velo = velo/sqrt(Tav/k/Tsoll)
      endif
c if it oscillates between bonded and not, reset more
      if(nfrag.eq.1) morestep=0

      if(nfrag.gt.1.and.tfrag.lt.1.d-6) tfrag=ttime/1000.

c check for fragmentation, EXIT section for frag runs
      if(it.gt.0)then
         call fragment(nuc,iat,pos,3.0d0,1,0,list)
         call fragmass(nuc,iat,list,mass,imass,nfrag,fragm,fragf,fragat)
         if(nfrag.gt.1)
     .   call intenergy(nuc,iat,list,mass,pos,velo,nfrag,fragT)
c probably an error
         if(nfrag.gt.6) goto 1000
c exit always immidiately if we have > nfragexit frags
         if(nfrag.gt.nfragexit) then
            write(*,8000)
     .      nstep,ttime,Epot,Ekin,Epot+Ekin,
     .      Eerror,nfrag,etemp,fragT(1:nfrag)
            write(*,9001)
            fragstate=1
            exit
         endif
         if(nfrag.eq.2)then
            fconst=fconst+1
         else
            fconst=0
         endif
c exit if nfrag=2 is constant for some time
         if(fconst.gt.1000) then
            write(*,8000)
     .      nstep,nstep*tstep/fstoau,Epot,Ekin,Epot+Ekin,
     .      Eerror,nfrag,etemp,fragT(1:nfrag)
            write(*,9002)
            fragstate=2
            exit
         endif
c add a few more cycles because fragmentation can
c directly proceed further and we don't want to miss this
         if(nfrag.ge.nfragexit) then
           morestep=morestep+1
           if(morestep.gt.more) then
              write(*,8000)
     .        nstep,ttime,Epot,Ekin,Epot+Ekin,
     .        Eerror,nfrag,etemp,fragT(1:nfrag)
              write(*,9003)
              fragstate=1
              exit
           endif
         endif
      endif

      enddo
c all done and nice
      mdok=.true.
      if(k.ge.nmax)then
         write(*,9004)
         fragstate=1
      endif

c error exit
1000  if( it.ge.0.and.it.lt.9999.and.(.not.restart) )close(1)
      if( it.eq.0) close(2)

c printed or used in main
      Tav =Tav/k
      Epav=Epav/k
      Ekav=Ekav/k
      aspin=avspin/kdump
      achrg=avchrg/kdump
      axyz =avxyz /kdump
      aTlast=aTlast/kdump

8000  format(i7,f8.0,F13.5,F9.4,F12.5,F8.4,1x,I1,F9.0,2x,5F7.0)
9001  format(/8x,'E X I T   M D  because of multiple fragments')
9002  format(/8x,'E X I T   M D  because nothing happens here anymore')
9003  format(/8x,'E X I T   M D  because of nfrag=nfragexit')
9004  format(/8x,'E X I T   M D  because of tmax has been reached')

      end

