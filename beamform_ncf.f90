!include 'fftw3.h'
!beamforming source origin analysis
!wangkai,2015/4/10
!junxie,2015/05/12
!       2016/02/12
!       2016/07/29 introduce sacio.mod
!parameter (nst=2500,nsmax=500000,nfmax=32768)
program beamforming
use sacio
implicit none
type(sac_head) :: sachead
integer  nst,nsmax,nfmax
parameter (nst=2500,nsmax=50000,nfmax=8192)
character(300) filein,par,dir,station(nst),stalst,outfile
integer nfreq,nsta,nframes,Ifmin,Ifmax,cor
integer npts,nlen,nerr,nsp
integer ferr,ista,i,j,npow
integer sublen,nsl,ifreq
integer narg,irr,is,nk
complex fseis(nfmax,15000)
complex sigout(nfmax)
real stlo(nst),stla(nst),freq(nfmax),sl(1000),theta(1000)
real er/6371.0/,pi/3.14159265/
real fmin,fmax,slmin,dsl
real sta/0/,sto/0/
real beg,delta,df,dt0
real sigin(nfmax)
real sig(nsmax)
logical ext
sig=0
sigout=cmplx(0,0)
narg=iargc()
if(narg.ne.1) then
  !write(*,*)'Usage: beamforming station.list param.dat cor_dir output'
   write(*,*)'Usage: beamforming param.dat'
   write(*,*)'param.dat is like:'
   write(*,*)'Station list'
   write(*,*)'NCF directory'
   write(*,*)'fmin,fmax slmin,dsl,nsl'
   stop
endif
! get filelist name
!-----------------
!call getarg (1,stalst)
call getarg (1,par)
!call getarg (3,dir)
!call getarg (2,outfile)
! get parameter from param.dat
!----------------------------
open(9,file=par)
read(9,'(a90)')stalst
read(9,'(a90)')dir
read(9,*)fmin,fmax,slmin,dsl,nsl
write(outfile,'(f5.3,"_",f5.3,".dat")')fmin,fmax
!write(*,*)outfile
! minimum frequency, maximum frequency, minimum slowness, delta_slowness, number of slowness
!----------------------------
close(9)
!write(*,*) fmin,fmax,slmin,dsl,nsl
open(10,file=stalst)
! read in station list
!----------------------------
nsta=0
do ista=1,nst
   read(10,*,iostat=ferr) station(ista),stlo(ista),stla(ista)
   if (ferr.lt.0)exit
   nsta=nsta+1              
   sta=sta+stla(ista)
   sto=sto+stlo(ista)
enddo
!write(*,*)'number of station is ',nsta
close(10)
      !nsta=nsta-1           ! number of stations
sta=sta/nsta
sto=sto/nsta
! get the new location of each stations taking the center as reference point
!----------------------------
do ista=1,nsta
   stlo(ista)=(stlo(ista)-sto)/180*pi*cos(sta/180*pi)*er
   stla(ista)=sin(stla(ista)*pi/180)/cos(sta/180*pi)*er
enddo
! read file with data and botain spectrum 
!----------------------------
irr=0
ista=0
do i=1, nsta-1
   do j=i+1,nsta
      ista=ista+1
      fseis(:,ista)=cmplx(0,0)
!     filein=trim(dir)//'/COR_'//trim(station(i))//'_'//trim(station(j))//'.SAC'
      write(filein,'(1a,"/",1a,"/l_COR_",1a,"_",1a,".SAC")')trim(dir),trim(station(i)),trim(station(i)),trim(station(j))
!                write(*,*)'read in file ',trim(filein)
      inquire(file=filein,exist=ext)
      if(.not.ext)cycle
      !call rsac1(trim(filein),sig,npts,beg,delta,nsmax,nerr)
      call read_sachead(filein,sachead,nerr)
      if(nerr.eq.-1)cycle
      call read_sac(filein,sig,sachead,nerr)
      if(nerr.eq.-1)cycle
      if(sachead%depmax.ne.sachead%depmax)cycle
      if(j.eq.2)then
         dt0=sachead%delta 
         npts=sachead%npts
         nsp=2
         npow=1
         do while (nsp.lt.npts)
            nsp=nsp*2
            npow=npow+1
         enddo
         nk=nsp/2+1
         df=1/dt0/nsp
         do ifreq=1,nsp/2+1
            freq(ifreq)=(ifreq-1)*df
         enddo
      endif
      if(sachead%delta.ne.dt0)cycle
      if(sachead%npts.ne.npts)cycle
      !npow=ceiling((log(dble(npts))/log(2.0)))
      !nsp=8192
      !npow=13
      !nsp=2**npow !determin power of FFT
      sigout=cmplx(0,0)
      sigout(1:npts)=cmplx(sig(1:npts),0)
      call clogc(npow,sigout,1,dt0)
      !call dfftw_plan_dft_r2c_1d(plan,nsp,sigin,sigout,FFTW_ESTIMATE)
      !call dfftw_execute_dft_r2c(plan,sigin,sigout)
      !call dfftw_destroy_plan(plan)
      fseis(1:nsp,ista)=sigout(1:nsp)
      !ista=ista+1
   enddo !  loop over station two
enddo    !  loop over station one
write(*,*)'hello'
write(*,*)'Total missing file is ',irr
!---------------------------------------------------
beg=0.0
!----------------------------------------------------
      !Nsta=ista-1            ! number of SAC traces
!      write(*,*)'number of sac traces is ',nsta
Nfreq=nsp/2+1;
Ifmin=int(fmin/df)+1   ! id of the minimum frequency
Ifmax=int(fmax/df)+1   ! id of the maximum frequency
do i=1,nsl
   sl(i)=slmin+(i-1)*dsl
!    write(*,*) sl(i)
enddo
do i=1,72
   theta(i)=(i-1)*5
!    write(*,*) theta(i)
enddo
! do beamforming
call fk(fseis,freq,stlo,stla,sl,theta,Nfreq,Nsta,Ifmin,Ifmax,nsl,outfile)
return
end

!========================================================
subroutine fk(fseis,freq,x,y,sl,theta,Nfreq,Nsta,Ifmin,Ifmax,nsl,outfile)
complex cj,tt
parameter(cj=(0.,1.))
real ::  pi=3.14159265
character(300) :: outfile
real, allocatable :: beam(:,:)
real x(1000),y(1000),sl(1000),theta(1000),beampower(72)
real projection(72,Nsta)
real isl,ifreq,omega
real freq(8192)
integer i,j,k,m,Ifmin,Ifmax
integer icc,iff,ith,itime,itc,Nfreq,Nsta,Nframes,Ntime,nsl
complex fseis(8192,15000)
complex rep(72,Nsta),temp1(72,Nsta),temp(72,72),repT(Nsta,72),cmat(Nsta,Nsta)
allocate(beam(72,Ifmax-Ifmin+1))
write(*,*) Nfreq,Nsta,Ifmin,Ifmax
do i=1,72                         ! loop over every angle
   do j=1,Nsta            ! loop over every stations
       projection(i,j)=x(j)*sin(theta(i)/180*pi)+y(j)*cos(theta(i)/180*pi);            
   enddo
enddo
write(*,*)'hello, see me? output file is ',outfile
!write(*,"(F15.7)") (projection(i,2),i=1,72)
open(100,file=trim(outfile))
do icc=1,nsl                      ! loop over every slowness
   isl=sl(icc)
   k=1
   do iff=Ifmin,Ifmax                 ! loop over frequency band
      ifreq=freq(iff)  
      omega=ifreq*2*pi
      m=1
      do i=1,Nsta
         do j=1,Nsta
            cmat(i,j)=(0.,0.)
            m=m+1
         enddo
      enddo
      m=1
      do i=1,Nsta-1
         do j=i+1,Nsta
            cmat(i,j)=fseis(iff,m)
            m=m+1
         enddo
      enddo
!         write(*,*) isl,ifreq,omega,theta(itime)
      rep=exp((cj*omega*isl)*projection) 
!         write(*,"(F15.7,F15.7)") (rep(1,i),i=1,Nsta)           
      repT=transpose(rep)
      temp1=matmul(rep,cmat)
      temp=matmul(temp1,conjg(rept))
        
      do ith=1,72
         beam(ith,k)=real(temp(ith,ith))**2
      enddo
      k=k+1
   enddo !end loop of iff
   beampower=10*log10(sum(beam,2))
!   write(*,"(F15.7)") (beam(i,2),i=1,72)
   do ith=1,72
      write(100,*) theta(ith),isl,beampower(ith)
   enddo 
enddo !end loop of icc
close(100)
deallocate(beam)
return
end 
