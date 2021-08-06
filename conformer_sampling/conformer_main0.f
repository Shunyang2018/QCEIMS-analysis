ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c important routine:
c author: Shunyang
c June 19th 2020
c get conformer using crest as start points for
c need set xtb and crest path before use
c https://stackoverflow.com/questions/44256136/read-a-file-with-an-unknown-number-rows-in-fortran
c https://www.tek-tips.com/viewthread.cfm?qid=461673
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program conformer
      implicit none
      include 'common1.inc'
      include 'ehtcommon.f'
      include 'setcommon.f'
      include 'omp_lib.h'
      real*8, dimension(:,:),allocatable :: xyz
      real*8, dimension(:),allocatable :: energy,prob
      integer, dimension(:),allocatable :: ntrajs, nmax0
      real*8 ,allocatable ::tmpxyz (:,:), velo(:,:)
      real*8 ,allocatable ::xyzr(:,:,:), velor(:,:,:)
      integer function omp_get_thread_num()
      character*80 atmp
      character*80 fout
      logical ex, ECP, metal3d, noecp, nometal,check
      character*2 tmpiat
      real*8 Temp, tmp, aTlast
      integer nuc, n, k, j, io, num, i, ntraj
      real*8  ams(107), Ekin
      integer,allocatable ::iat (:)
      real*8 ,allocatable ::mass(:),imass(:)
      real*8 ,allocatable ::eps  (:,:)
      real*8 ,allocatable ::mopop(:,:)
      integer mchrg,na,nb,prog,iseed(1),iprog,nfragexit,maxsec,edistri
      real*8 ehomo, eimp0,exc,ieeel,tstep,tmax,etempin
      real*8 Tinit,iee_a,iee_b,eimpw,fimp,trelax,hacc,btf,ieeatm
      real*8 lowerbound,upperbound,Iacc, dum, fstoau,pmax,dums,Tav,edum
      real*8 tta,snorm,x,tsoll,tadd
      integer method,scan,unity,nrun,mo1,mo2,m
      real*8  gaus(1000)
      real*8 ,allocatable ::modum(:),velof(:)
c  atomic masses
      data  ams /  1.00790d0,  4.00260d0,  6.94000d0,  9.01218d0,
     110.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0,
     220.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0,
     330.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0,
     440.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0,
     554.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0,
     665.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0,
     779.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0,
     891.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0,
     9102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0,
     1118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0,
     2132.9054d0, 137.3300d0, 15*0.000d0, 178.4900d0, 180.9479d0,
     3183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0,
     4196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0,
     518*0.000d0,   0.0000d0,  5*0.000d0/


      aTlast=3000.0
c timesteps in au
      tstep=tstep*fstoau
c the "70" eV
      eimp0=eimp0/27.2113957

      call input(prog,tstep,tmax,ntraj,iseed(1),etempin,Tinit,
     .           iee_a,iee_b,eimp0,eimpw,fimp,iprog,
     .           trelax,hacc,nfragexit,maxsec,edistri,btf,ieeatm,
     .           method,scan,lowerbound,upperbound,metal3d,Iacc,
     .           ECP,unity,noecp,nometal,gfnver)
c     xtb optimization
      fout = 'xtb.log'
      write(atmp,'(''xtb coord --opt vtight > ''a)') trim(fout)
      !call system(atmp)
      inquire(file='xtbopt.tmol',exist=ex)
      if(.not.ex) stop 'fatal QC error. Must stop!'
      if(.not.ex) stop 'coordinates optimization failed'

c       crest sampling
      fout = 'crest.log'
      write(atmp,'(''crest xtbopt.tmol -temp ''F8.1'' > ''a)') Temp,
     .  trim(fout)
      !call system(atmp)
      inquire(file='crest_conformers.xyz',exist=ex)
      if(.not.ex)then
        stop 'crest calculation failed'
      else
        write(*,*)'getting conformers from crest_conformers.xyz...'
        open(unit=2,file='crest_conformers.xyz')
        n = 0
        read(2,*) nuc
        n=n+1
c get full length of the xyz file
        do
          read(2,*,iostat=io)
          if (io/=0) EXIT
          n = n+1
c          write(*,*)'n', n
        end do
        num = n/(nuc+2) ! num is number of conformers
        write(*,'(a, i, a)') 'have read', num, '  conformers'


        rewind(2)! important!!! read from very begining
c read xyz and energy
      if (num .gt. 100) num=100
      write(*,*)'maximum',num,'has been read!'
      allocate(xyz(3,num*nuc),energy(num),prob(num),ntrajs(num),
     .nmax0(num))
      do j = 1, num

        read(2,*,iostat=io) nuc ! important! must know whether end of file at this step,
                              ! because the do loop according to nuc won't be excess
c        if (io/=0) Exit !read to end of file
        read(2,*) tmp
        energy(j) = tmp ! reading energy for boltzmann distribution

        do k=1,nuc
          read(2, *,iostat=io) tmpiat, (xyz(i,(j-1)*nuc+k),i=1,3) !nuc not num!!!!
c          write(*,*) k,j, (xyz(i,(j-1)*num+k),i=1,3)
c          write(*,*)io
        end do

      end do

        endif
      write(*,*)'reading conformers done!'
        close(2)


c calculating boltzmann distribution
      write(*,*)'calculating boltzmann distribution'

      write(*,*)'energy',energy
      call boltz(2,num,Tinit,energy-energy(1),prob)

      write(*,*)'Temperature', Tinit


c set trajectories
      ntraj = nuc *25
      ntrajs = int(prob*ntraj)
      ntraj = sum(ntrajs)
      write(*,*) 'total ntraj', ntraj
      write(*,*) 'ntraj of different conformers', ntrajs
      write(*,*)'probability of each conformer',prob
      nmax0 = ntrajs*50
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      do k = 1,num
        write(*,*)'md equilibrate',num
        i=1
        write(*,*)'md 2',OMP_GET_THREAD_NUM()
      end do
!$OMP END DO
!$OMP END PARALLEL

      end program
