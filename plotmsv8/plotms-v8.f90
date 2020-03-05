! ============================================================================
! Name        : plotms.f90
! Author      : S. Grimme, modifications T. Kind (FiehnLab 2013)
! Version     : 2.4  (Feb 19 2014)
! Copyright   :
! Description : plotms from QCEIMS
! ============================================================================

!====================================================================
!     Original PlotMS Program for extraction of mass spectra
!     from quantum mechanical simulations used within QCEIMS
!
!                      *********************************************
!                      *                S. Grimme                  *
!                      * Mulliken Center for Theoretical Chemistry *
!                      *             Universitaet Bonn             *
!                      *                  2008-13                  *
!                      *********************************************
!
!     Please cite as:
!     S.Grimme, Angew.Chem.Int.Ed. 52 (2013) 6306-6312.
!
!     Adapted version:
!     Tobias Kind - FiehnLab 2013
!
!     Changes:
!     1) Output MSP or JCAMP formated unit mass mass spectra (requires exp.dat and mass.agr)
!     2) Elements Br and W added
!     3) Chlorine masses changed to accurate masses
!     4) Output of exact masses
!
!     Adapted version:
!     Shunyang Wang - FiehnLab 2020
!
!
!====================================================================




      program plotms
      use dictionary_m ! get hash table
!     treat i,j,k,l,m,n as integers and all other cariables as real
      implicit none
!    declare parameters in dictionary
      type(dictionary_t) :: dict, tdict
      character(len=16), allocatable :: key_list(:), tkey_list(:)
      integer dictsize,s
      character*16 charstrk, charstrv
      real*16 ttmp, tmp, floatv

      integer n,i,j,k,kk,kkk,kkkk,nn,natot,nsig
      integer nrnd,ndim,imax,imin,nagrfile
! maximum mass = 10000
      real*8  iexp (2,10000)
      integer iat  (10000)
      integer nat  (10000)
      integer idum (10000)
      integer mass (1000)
      integer isec,jsec,ial,jal(0:10),nf,irun
      real*8  mint (1000)
! TK  tmass contains abundances (real) for each unit mass (i)
      real*8  tmass(10000)
      real*8  checksum2(50000)
      real*4, allocatable :: rnd(:,:)
! TK  cthr and cthr contain ions (counts) with different charges
!     tmax the maximum abundance factor over the whole spectrum
      real*8 xx(100),tmax,r,rms,norm,chrg,cthr,cthr2
      real*8 chrg1,chrg2,dum,checksum,cw
      real   snorm
      logical ex,sel,echo,exdat,mpop

! TK  fname=<qceims.res> or result file, xname contains the mass.agr plot file
      character*80 arg(10),line,fname,xname
      character*255 atmp

      ! TK number of peaks in in-silico spectra, needed for JCAMP-DX export
      integer numspec
      ! initialize the hash TABLE
      dictsize=1024
      call dict%init(dictsize)

      iexp =0
      mpop =.false.
      echo =.false.
      isec =0
      cthr =1.d-3
      cw   =0
      norm =1.0
      cthr2=1.01
      nagrfile=410

! edit this path name to some standard xmgrace plot file
! TK changed to direct working path

      xname='mass.agr'
      fname=''

      ! TK start loop reading arguments
      do i=1,9
         arg(i)=' '
         call getarg(i,arg(i))
      enddo
      ! TK end loop reading arguments


      ! TK start mega loop processing aruments
      do i=1,9

! TK comand line parameters
! TK -a no idea (related to ion charge count)
! TK -p print spectra "mass % intensity  counts   Int. exptl" to stdout;
!       with "Int. exptl" (experimental) taken from exp.dat but not all exp peaks are exported
!       if no theoretical counterpart exists
! TK -f filename or  -f <name_of_res_file>
! TK -t no idea
! TK -w no idea
! TK -s no idea

         if(index(arg(i),'-a').ne.0)  cthr=-1000.
         if(index(arg(i),'-p').ne.0)  echo=.true.
         if(index(arg(i),'-f').ne.0)  fname=arg(i+1)
         if(index(arg(i),'-t').ne.0) then

            ! TK call subroutine READL
            call readl(arg(i+1),xx,nn)
            cthr=xx(1)
         endif
         if(index(arg(i),'-w').ne.0) then
            call readl(arg(i+1),xx,nn)
            cw=xx(1)
         endif
         if(index(arg(i),'-s').ne.0) then
            call readl(arg(i+1),xx,nn)
            isec=int(xx(1))
         endif
      ! TK end mega loop processing aruments
      enddo

! fname contains the results from each calculation or the temporary result tmpqceims.res
! xname contains the xmgrace plot file

      if(fname.eq.'')then
        fname='qceims.res'
        inquire(file=fname,exist=ex)
        if(.not.ex) fname='tmpqceims.res'
      endif

      inquire(file=fname,exist=ex)
      if(.not.ex) stop 'res file does not exist'

      write(*,*) 'QCEIMS output reader PLOTMS'
      write(*,*) 'V 2.2, Nov 2013'
      write(*,*)
      write(*,*) 'xmgrace file body ',trim(xname)
      write(*,*) 'Reading ... ', trim(fname)

      if(cthr.ge.0)then
      write(*,'( &
      '' couting ions with charge from '',f4.1,'' to '',f4.1)')   cthr,cthr2
      else
      write(*,'( &
      '' counting all fragments with unity charge (frag. overview)'')')
      endif
      if(cw.gt.0)then
      write(*,'( &
     '' broadening the charges by an SD, wdth :'',f4.1)')cw*100.
      endif

      if(isec.ne.0)then
      write(*,*) &
      'Taking only secondary, tertiary ... fragmentations ',isec
      endif

      tmass = 0
      call tdict%init(dictsize)
      call random_seed()

! read file once
      i=1
      ndim=0

      ! TK contains qceims.res as standard option
      ! TK example output of qceims.res with 13 variables (in this case)
      ! 0.1986979   21 1 1    4     1  5   6  6   7  5   8  1

      !  other example from qceims.out
      !  trajectory          100           1
      !  mass                 formula              q pop   spin    q IPB  diss time (ps)
      !  M=121.12               H5C7O2   100   1   0.361   0.985   0.000    0.426
      !  M=73.19               H9C3SI1   100   1   0.639   0.015   1.000    0.426 ~

      !  is coded in qceims.res as
      !             trj # sub  #elem H  5   C  7   O  2
      !  0.0000088  100 1 1    3     1  5   6  7   8  2

      !             trj # sub  #elem H  9   C  3  Si  1
      !  0.9999912  100 2 1    3     1  9   6  3  14  1


      ! TK irun is optimized out during compile time
      ! example qceims.res
      ! 0.0000000   13 1 1    2     6  1   8  2

      ! chrg2 = 0
      ! irun = 13
      ! jsec = 1
      ! nf = 1
      ! k = 2 // kk = 3 (?)
      ! iat(kk) = 0
      ! nat(kk) = 0
      ! k = 2

      open(unit=1,file=fname)
 10   read(1,*,end=100)chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)

      ! TK just for debug purposes
      write(*,*) chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
      ! TK end debug purposes

      ! TK normal program
      if(isec.gt.0.and.isec.ne.jsec) goto 10
      if(chrg2.gt.cthr) then
         natot=sum(nat(1:k))
         if(natot.gt.ndim)ndim=natot
      endif
      i=i+1
      goto 10


100   continue
      n=i-1
      ! TK n contains number of fragments and ndim contains number of atoms(?) max
      write(*,*) n,' fragments with ',ndim,' atoms max.'
      close(1)

! initialize the random number array (efficiency)
      nrnd=50000
      allocate(rnd(ndim,nrnd))
      do i=1,nrnd
         do j=1,ndim
            call random_number(r)
            rnd(j,i)=r
         enddo
      enddo

! read it again
      write(*,*) 'Computing ...'

! TK contains qceims.res as standard option
      open(unit=1,file=fname)
      imin=100000
      i=1
      ial=0
      jal=0
      checksum =0
      checksum2=0
 11   read(1,*,end=101)chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
      sel=.false.
      chrg=chrg2
      checksum2(irun)=checksum2(irun)+chrg2
      if(cthr.lt.0)chrg=1.0
      if(chrg.gt.cthr)then
         sel=.true.
         ial=ial+1
         jal(jsec)=jal(jsec)+1
      endif
      if(isec.gt.0.and.isec.ne.jsec) sel=.false.
      if(sel)then
        natot=sum(nat(1:k))
! all types of atoms
        kkkk=0
        do kk=1,k
! all atoms of this type
          do kkk=1,nat(kk)
             kkkk=kkkk+1
             idum(kkkk)=iat(kk)
          enddo
        enddo
        checksum=checksum+chrg

! compute pattern, nsig signals at masses mass with int mint
        write(*,*)'calculating traj',irun,'step',jsec,'frag',nf
        call isotope(natot,idum,ndim,nrnd,rnd,mass,mint,nsig,dict,key_list)
        if(cw.gt.1.d-6)chrg=chrg+cw*chrg*snorm()
        call dict%show(key_list)
        s = SIZE(key_list)
! SW debug print
      !  write(*,*)'check peaks'
      !  do i = 1, s
      !    write(*,*) key_list(i)
      !  end do
      !  if (s .ne. nsig) write(*,*)s,'wrong with round up',nsig
        do k=1,nsig
         charstrk=key_list(k)
         tmass(mass(k))=tmass(mass(k))+mint(k)*chrg

! SW accurate mass
        !write(*,*)'key value used in floatv function',charstrk
         tmp = floatv(charstrk,dict)

         if (tdict%get(charstrk) .eq. ' ') THEN
           write(charstrv,'(F16.6)')tmp
           !write(*,*)'charstrv',charstrv
           call tdict%set(charstrk, charstrv)
         else
           ttmp = floatv(charstrk,tdict) ! total peak intensity
          ! write(*,*)'ttmp*chrg',ttmp
           ttmp = ttmp +tmp
           WRITE (charstrv,'(F16.6)')ttmp
           call tdict%set(charstrk, charstrv)
         endif
        enddo
        deallocate(key_list)
        i=i+1
      endif
      goto 11
101   continue
      write(*,*) n,' (charged) fragments done.'
      write(*,*)
      write(*,*) 'checksum of charge :',checksum
      write(*,*)

      k=0
      do i=1,50000
         if(abs(checksum2(i)).gt.1.d-6)k=k+1

         if(abs(checksum2(i)).gt.1.d-6.and. &
            abs(checksum2(i)-1.).gt.1.d-3)then
            write(*,*) 'checksum error for trj', i,' chrg=',checksum2(i)
         endif
      enddo
      write(*,*) k,' successfull runs.'

      write(*,*)
      write(*,*) 'ion sources:'
      write(*,'(''% inital run'',f6.1)') &
                100.*float(jal(1))/float(ial)
      do j=2,8
      write(*,'(''%        '',i1,''th'',f6.1)')j, &
                100.*float(jal(j))/float(ial)
      enddo
      write(*,*)

! read exp.
      imax=0
      imin=0
      inquire(file='exp.dat',exist=exdat)
      if(exdat) then
         write(*,*) 'Reading exp.dat ...'
         open(unit=4,file='exp.dat')

! TK (a) is edit descriptor for character strings
20       read(4,'(a)',end=200)line
         if(index(line,'##MAXY=').ne.0)then
            line(7:7)=' '
            call readl(line,xx,nn)
            norm=xx(nn)/100.0d0
         endif
         if(index(line,'##PEAK TABLE').ne.0)then
            kk=0
30          read(4,'(a)',end=200)line

            ! TK JCAMP DX for MS data has "##END="
            if(index(line,'##END').ne.0)goto 200
            do k=1,80
               if(line(k:k).eq.',')line(k:k)=' '
            enddo
            call readl(line,xx,nn)
            do k=1,nn/2
               kk=kk+1
               iexp(1,kk)=xx(2*k-1)
               iexp(2,kk)=xx(2*k)
               if(iexp(1,kk).gt.imax) imax=iexp(1,kk)
               if(iexp(1,kk).lt.imin) imin=iexp(1,kk)
            enddo
            goto 30
         endif
         goto 20
200      continue
         close(4)
      endif

      imin=max(imin,10)
      j=0
      do i=10,10000
         if(tmass(i).ne.0)j=i
      enddo
      imax=max(j,imax)

      tmax=maxval(tmass(10:10000))
      ! TK number of peaks computed in theoretical spectrum
      ! idint(tmax) is not related to the number of spectral peaks
      write(*,*)'Theoretical counts in 100 % signal:',idint(tmax)

      if(echo)then
      write(*,*)'mass % intensity  counts   Int. exptl'
      do i=10,10000
         if(tmass(i).ne.0)then
            dum=0
            do j=1,kk
               if(int(iexp(1,j)).eq.i)dum=iexp(2,j)
            enddo
            write(*,'(i8,F8.2,i8,F8.2)') &
            i,100.*tmass(i)/tmax,idint(tmass(i)),dum/norm
         endif
      enddo
      endif

! when the template exists
!TK Fortran runtime error: File already opened in another unit
!TK At line 285 of file src/plotms.f90 (unit = 3, file = '')

      inquire(file=xname,exist=ex)
      if(ex)then
      open(unit=2,file=xname)

! TK original mass.agr a preconfigured file
      open(unit=7,file='mass-result.agr')

! my xmgrace file mass.agr has nagrfile lines
      write(*,*) 'Writing mass-result.agr ...'
      do i=1,nagrfile
         read (2,'(a)')line
         if(index(line,'world 10').ne.0)then
            write(line,'(''@    world'',i3,'', -105,'' &
                           ,i3,'', 100'')')imin,imax+5
         endif
         if(index(line,'xaxis  tick major').ne.0)then
            line='@    xaxis  tick major 20'
            if(imax.gt.200) line='@    xaxis  tick major 50'
         endif
         if(index(line,'@    s1 symbol size').ne.0)then
            line='@    s1 symbol size 0.160000'
            if(imax.gt.200) line='@    s1 symbol size 0.100000'
         endif
         if(index(line,'@    s2 symbol size').ne.0)then
            line='@    s2 symbol size 0.160000'
            if(imax.gt.200) line='@    s2 symbol size 0.100000'
         endif
         write(7,'(a)')line
      enddo

! write only masses > 10
! TK tmass contains theoretical masses, unit 7 = mass-result.agr

      do i=10,10000
         if(tmass(i).ne.0)then
            write(7,*) i,100.*tmass(i)/tmax
         endif
      enddo
      write(7,*)'&'

! TK establish again if exdat (experimental JCAMPDX) exists
! TK @target G0.S1 means upper theory spectrum in current implementation
! TK @target G0.S2 means lower experimental spectrum
! TK iexp(1,i) - contains the exp masses and  iexp(2,i) contains the abundance
! TK kk contains number of experimental spectra

      if(exdat)then
      write(7,*)'@target G0.S2'
      write(7,*)'@type bar'
      do i=1,kk
         write(7,*) iexp(1,i),-iexp(2,i)/norm
      enddo
      write(7,*)'&'
      endif
      endif

! TK here we export the spectrum to JCAMP-DX, subroutine would be better for
! TK allowing later code merges or inclusion of new and missing elements such as
! TK currently unit mass only
    ! SW also generate accurate mass
    ! TK open file as JCAMP-DX MS result file, open new or replace
    write(*,*)'open file 11'
    open(unit=11, file='result.jdx', STATUS="REPLACE")
    write(*,*)'open file 111'
    open(unit=111, file='accuratemass.jdx',STATUS="REPLACE")
    ! TK minimum JCAMP-DX header
    write(11,"(A)")'##TITLE=Theoretical in-silico spectrum (QCEIMS)'
    write(11,"(A)")'##JCAMP-DX=4.24'
    write(11,"(A)")'##DATA TYPE=MASS SPECTRUM'
! SW
    write(111,"(A)")'##TITLE=Theoretical in-silico spectrum (QCEIMS)'
    write(111,"(A)")'##JCAMP-DX=4.24'
    write(111,"(A)")'##DATA TYPE=MASS SPECTRUM'
    !write(11,*)'##MOLFORM=C13 H22 O3 Si2'
    !write(11,*)'##MW=282'

    write(11,"(A)")'##XUNITS=M/Z'
    write(11,"(A)")'##YUNITS=RELATIVE INTENSITY'
!SW
    write(111,"(A)")'##XUNITS=M/Z'
    write(111,"(A)")'##YUNITS=RELATIVE INTENSITY'
    !write(11,*)'##XFACTOR=1'
    !write(11,*)'##YFACTOR=1'
    !write(11,*)'##FIRSTX=15'
    !write(11,*)'##LASTX=281'
    !write(11,*)'##FIRSTY=20'
    !write(11,*)'##MAXX=281'
    !write(11,*)'##MINX=15'
    !write(11,*)'##MAXY=9999'
    !write(11,*)'##MINY=10'

    ! TK calculate number of in-silico spectra
    ! TK new: numspec = idint(tmax) ** not related to number ofr spectral peaks
     numspec = 0
     do i=10,10000
         if(tmass(i).ne.0)then
            numspec = numspec + 1
         endif
     enddo


    write(11,'(A, I0)')'##NPOINTS=' ,numspec
    call tdict%show(tkey_list)
    s = SIZE(tkey_list)
    write(111,'(A, I0)')'##NPOINTS=' ,s

    ! TK ##XYDATA=(XY..XY) 1 designates one line per m/z abd pair
    ! TK separated by comma
    ! TK 1.500000e+001 ,5.000000e+000

    write(11,"(A)")'##PEAK TABLE=(XY..XY) 1'
    write(111,"(A)")'##PEAK TABLE=(XY..XY) 1'
    !15,20 26,10 27,20 29,50

    write(*,*)'write result.jdx....'
    do i=10,10000
         if(tmass(i).ne.0)then
           ! unit mass to exact mass
           ! write(11,"(I4, I8)") i, int(100.*tmass(i)/tmax)
           write(11,"(F12.6, I8)") real(i), int(100.*tmass(i)/tmax)
         endif
    enddo
    write(*,*)'write accurate mass...'
    do i = 1, s
      write(111,*) tkey_list(i), tdict%get(tkey_list(i))
    end do

    write(11,"(A)")'##END='
    write(111,"(A)")'##END='

    ! TK close the JCAMP-DX file
    close(11)
    close(111)

! compute deviation exp-theor.
! TK here we potentially have to use Mass Spec related terms, such as dot product
! From Stein and Scott

      if(exdat)then
      rms=0
      nn=0
      kkk=0
      kkkk=0
      do i=1,kk
         k=iexp(1,i)
         if(iexp(2,i)/norm.gt.5.0)then
            r=100.*tmass(k)/tmax-iexp(2,i)/norm
            kkk=kkk+1
            if(100.*tmass(k)/tmax.gt.2.5)kkkk=kkkk+1
            rms=rms+abs(r)
         endif
      enddo
      write(*,*)'MAD(exptl./theor.) = ',rms/kkk
      write(*,*)'# exptl. > 5 %     = ',kkk
      write(*,*)'% correctly found  = ',100*float(kkkk)/float(kkk)
      endif

      close(1)
      close(2)
      close(3)
      close(7)
      end
! compute the isotopic mass distribution for
! given chemical formula (nat atoms with OZ it())
! returns nsig masses in mass() with probability
! mint()
! this is a quick and dirty MC algorithm and its run-time depends
! critically on the number of trials nrnd
! TK currently 10 elements covered:
! TK 1H, 6C, 7N, 8O, 9F, 14Si, 15N, 16S, 17Cl, 26Fe
! TK increased mass acuracy values, currently unit masses
! TK See: ATOMIC WEIGHTS OF THE ELEMENTS: REVIEW 2000 (IUPAC Technical Report)
! TK Pure Appl. Chem., Vol. 75, No. 6, pp. 683�800, 2003.
! SW calculate accurate mass
      subroutine isotope(nat,it,ndim,nrnd,rnd,mass,mint,nsig,dict,key_list)
      use dictionary_m
      implicit none
      character(len=16), allocatable :: key_list(:)
      character(len=16) charstrk, charstrv
      type(dictionary_t):: dict
      real*16 intensity
      integer nat,it(*),nsig,nrnd,ndim, s, dictsize
      integer mass(*)
      real*8  mint(*)
      real*4  rnd(ndim,nrnd)

      real*4  prob(200,4),massiso(200,4),p1,p2,x,xmass,examass
      integer niso(200)
      integer nmass(10000)
      integer n,i,j,iso,imass,isum,k,iti
      real*8 r,xm

      niso=0
      prob=0
      massiso=0
! SW accuatemass
 ! TK 1H
      niso(1)=2
      prob(1,1)=0.0115
      prob(1,2)=99.9885
      massiso(1,1)=2.014102
      massiso(1,2)=1.007825

 ! TK 6C
      niso(6)=2
      prob(6,1)=1.07
      prob(6,2)=98.93
      massiso(6,1)=13.003355
      massiso(6,2)=12.000000

 ! TK 7N
      niso(7)=2
      prob(7,1)=0.368
      prob(7,2)=99.632
      massiso(7,1)=15.000109
      massiso(7,2)=14.003074

 ! TK 8O (Oxygen)
      niso(8)=3
      prob(8,1)=0.038
      prob(8,2)=0.205
      prob(8,3)=99.757
      massiso(8,1)=16.999132
      massiso(8,2)=17.999160
      massiso(8,3)=15.994915

 ! TK 9F
      niso(9)=1
      prob(9,1)=100.
      massiso(9,1)=18.998403

 ! TK 14Si
      niso(14)=3
      prob(14,1)=92.223
      prob(14,2)=4.685
      prob(14,3)=3.092
      massiso(14,1)=27.976927
      massiso(14,2)=28.976495
      massiso(14,3)=29.973770

! SW 15P
 ! TK 15N
      niso(15)=1
      prob(15,1)=100.
      massiso(15,1)=30.973762

 ! TK 16S
      niso(16)=4
      prob(16,1)=0.02
      prob(16,2)=0.76
      prob(16,3)=4.29
      prob(16,4)=94.93
      massiso(16,1)=35.967081
      massiso(16,2)=32.971458
      massiso(16,3)=33.967867
      massiso(16,4)=31.972071

 ! TK 17Cl - changed to isotopic masses and abundances OK
      niso(17)=2
      prob(17,1)=75.76
      prob(17,2)=24.24
      massiso(17,1)=34.968853
      massiso(17,2)=36.965903

! TK 26Fe
      niso(26)=4
      prob(26,1)=5.845
      prob(26,2)=91.754
      prob(26,3)=2.119
      prob(26,4)=0.2
      massiso(26,1)=53.939
      massiso(26,2)=55.934
      massiso(26,3)=56.935
      massiso(26,4)=57.933

! TK added 35Br
      niso(35)=2
      prob(35,1)=50.69
      prob(35,2)=49.31
      massiso(35,1)=78.91833
      massiso(35,2)=80.91629


      prob = prob * 0.01
      dictsize=1024
      call dict%init(dictsize)
! TK mass currently only loops to element 36  (Krypton)
      do i=1,36
         xm=0
         do j=1,niso(i)
            xm=xm+prob(i,j)
         enddo
         if(niso(i).gt.0.and.abs(xm-1.).gt.0.01) &
         stop 'internal isotope error 1'
      enddo

      do i=1,nat
         if(it(i).gt.100)then
            niso   (it(i))  =1
            prob   (it(i),1)=1.
            massiso(it(i),1)=it(i)-100.
         endif
      enddo

      do i=1,nat
         if(niso(it(i)).eq.0) stop 'internal isotope error 2'
      enddo

      nmass=0
!      nrnd=200
      do n=1,nrnd
      xmass=0
      do i=1,nat
         iti=it(i)
         r=rnd(i,n)
         p1=0.0
         p2=prob(iti,1)
         do iso=1,niso(iti)
            if(r.ge.p1.and.r.le.p2)then
               x=massiso(iti,iso)
               exit
            endif
            p1=p2
            p2=p2+prob(iti,iso+1)
         enddo
         xmass=xmass+x
         !write(*,*)'exact mass',xmass
      enddo
      !SW mass and intensity
      imass=int(xmass)
      nmass(imass)=nmass(imass)+1
      !SW exact mass and intensity
      write(charstrk,'(F16.6)')xmass
    !  write(*,*) charstrk,'<key vaule>',dict%get(charstrk)
      !WRITE (STRING,*) variable.
      if (dict%get(charstrk) .eq. ' ') THEN
        charstrv='0.00002'
        call dict%set(charstrk, charstrv)
      else
        charstrv = dict%get(charstrk)
        read (charstrv,*)intensity
        !write (*,*)'intensity',intensity
        intensity = (intensity+0.00002) !normalization
        write(charstrv,'(F16.6)')intensity
      !  write (*,*)'charstrv',charstrv
        call dict%set(charstrk, charstrv)
      endif
      enddo

      isum=sum(nmass)
      k=0
      do i=1,2000
         if(nmass(i).gt.0)then
            k=k+1
            mass(k)=i
            mint(k)=float(nmass(i))/float(isum)

         endif
      enddo


      nsig=k
      end


! TK sub function as explained below
      REAL FUNCTION snorm()

!C**********************************************************************C
!C                                                                      C
!C                                                                      C
!C     (STANDARD-)  N O R M A L  DISTRIBUTION                           C
!C                                                                      C
!C                                                                      C
!C**********************************************************************C
!C**********************************************************************C
!C                                                                      C
!C     FOR DETAILS SEE:                                                 C
!C                                                                      C
!C               AHRENS, J.H. AND DIETER, U.                            C
!C               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             C
!C               SAMPLING FROM THE NORMAL DISTRIBUTION.                 C
!C               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          C
!C                                                                      C
!C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  C
!C     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  C
!C                                                                      C
!C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
!C     SUNIF.  The argument IR thus goes away.                          C
!C                                                                      C
!C**********************************************************************C
!C
      DIMENSION a(32),d(31),t(31),h(31)
!C
!C     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
!C     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
!C

      DATA a/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,.1970991, &
           .2372021,.2776904,.3186394,.3601299,.4022501,.4450965, &
           .4887764,.5334097,.5791322,.6260990,.6744898,.7245144, &
           .7764218,.8305109,.8871466,.9467818,1.009990,1.077516, &
           1.150349,1.229859,1.318011,1.417797,1.534121,1.675940, &
           1.862732,2.153875/
      DATA d/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243, &
           .1899108,.1812252,.1736014,.1668419,.1607967,.1553497, &
           .1504094,.1459026,.1417700,.1379632,.1344418,.1311722, &
           .1281260,.1252791,.1226109,.1201036,.1177417,.1155119, &
           .1134023,.1114027,.1095039/
      DATA t/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2, &
           .7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1, &
           .1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1, &
           .2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1, &
           .4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1, &
           .9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980, &
           .5847031/
      DATA h/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1, &
           .4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1, &
           .4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1, &
           .4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1, &
           .5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1, &
           .8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016, &
           .7010474/

!C
   10 call random_number(u)
      s = 0.0
      IF (u.GT.0.5) s = 1.0
      u = u + u - s
   20 u = 32.0*u
      i = int(u)
      IF (i.EQ.32) i = 31
      IF (i.EQ.0) GO TO 100
!C
!C                                START CENTER
!C
   30 ustar = u - float(i)
      aa = a(i)
   40 IF (ustar.LE.t(i)) GO TO 60
      w = (ustar-t(i))*h(i)
!C
!C                                EXIT   (BOTH CASES)
!C
   50 y = aa + w
      snorm = y
      IF (s.EQ.1.0) snorm = -y
      RETURN
!C
!C                                CENTER CONTINUED
!C
   60 call random_number(u)
      w = u* (a(i+1)-aa)
      tt = (0.5*w+aa)*w
      GO TO 80

   70 tt = u
      call random_number(ustar)
   80 IF (ustar.GT.tt) GO TO 50
   90 call random_number(u)
      IF (ustar.GE.u) GO TO 70
      call random_number(ustar)
      GO TO 40
!C
!C                                START TAIL
!C
  100 i = 6
      aa = a(32)
      GO TO 120

  110 aa = aa + d(i)
      i = i + 1
  120 u = u + u
      IF (u.LT.1.0) GO TO 110
  130 u = u - 1.0
  140 w = u*d(i)
      tt = (0.5*w+aa)*w
      GO TO 160

  150 tt = u
  160 call random_number(ustar)
      IF (ustar.GT.tt) GO TO 50
  170 call random_number(u)
      IF (ustar.GE.u) GO TO 150
      call random_number(u)
      GO TO 140

      END

!     *****************************************************************
! TK  This routine handles
! TK  Called by: Main program
! TK  USES:      READAA

      SUBROUTINE READL(A1,X,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) A1
      DIMENSION X(*)
      I=0
      IS=1
  10  I=I+1
      X(I)=READAA(A1,IS,IB,IE)
      IF(IB.GT.0 .AND. IE.GT.0) THEN
                                IS=IE
                                GOTO 10
      ENDIF
      N=I-1
      RETURN
      END


! *******
! SW function to get value as float
      FUNCTION floatv(charstrk, d)
      use dictionary_m
      implicit NONE
      CHARACTER*16 charstrk, charstrv
      Real*16 floatv,intensity
      type(dictionary_t) :: d
      character(len=16), allocatable :: key_list(:)
      charstrv = d%get(charstrk)
!      write(*,*)'charstrk',charstrk
!      write(*,*)'charstrv',charstrv
      read (charstrv,*) intensity
      floatv = intensity
      return
    end function floatv

!     *****************************************************************
! TK  This function handles
! TK  Called by: Main program
! TK  USES:      no other sub

      FUNCTION READAA(A,ISTART,IEND,IEND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 READAA
      CHARACTER*(*) A
      NINE=ICHAR('9')
      IZERO=ICHAR('0')
      MINUS=ICHAR('-')
      IDOT=ICHAR('.')
      ND=ICHAR('D')
      NE=ICHAR('E')
      IBL=ICHAR(' ')
      IEND=0
      IEND2=0
      IDIG=0
      C1=0
      C2=0
      ONE=1.D0
      X = 1.D0
      NL=LEN(A)
      DO 10 J=ISTART,NL-1
         N=ICHAR(A(J:J))
         M=ICHAR(A(J+1:J+1))
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO .OR. M.EQ.IDOT)) GOTO 20
   10 CONTINUE
      READAA=0.D0
      RETURN
   20 CONTINUE
      IEND=J
      DO 30 I=J,NL
         N=ICHAR(A(I:I))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C1=C1*10+N-IZERO
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN
            ONE=-1.D0
         ELSEIF(N.EQ.IDOT) THEN
            GOTO 40
         ELSE
            GOTO 60
         ENDIF
   30 CONTINUE
   40 CONTINUE
      IDIG=0
      DO 50 II=I+1,NL
         N=ICHAR(A(II:II))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C2=C2*10+N-IZERO
            X = X /10
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN
            X=-X
         ELSE
            GOTO 60
         ENDIF
   50 CONTINUE

!C
!C PUT THE PIECES TOGETHER
!C

   60 CONTINUE
      READAA= ONE * ( C1 + C2 * X)
      DO 55 J=IEND,NL
         N=ICHAR(A(J:J))
         IEND2=J
         IF(N.EQ.IBL)RETURN
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
      RETURN

   57 C1=0.0D0
      ONE=1.0D0
      DO 31 I=J+1,NL
         N=ICHAR(A(I:I))
         IEND2=I
         IF(N.EQ.IBL)GOTO 70
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
         IF(N.EQ.MINUS)ONE=-1.0D0
   31 CONTINUE
   61 CONTINUE
   70 READAA=READAA*10**(ONE*C1)
      RETURN
      END
