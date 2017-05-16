         program beamforming
 include 'fftw3.f'
! wangkai,2015/4/10
 parameter (nst=2000,nsmax=500000,nfmax=32768)

 character*300 filein 
 real sig(nsmax)
 real beg,delta,df
 integer npts,nlen,nerr,nsp
 integer  ferr,ista,i,j,npow
 
 integer sublen,nsl
 real fmin,fmax,slmin,dsl
 complex fseis(nfmax,nst)
 real*8 sigin(nfmax)
 double complex sigout(nfmax)
 integer*8 plan
  
 real stlo(nst),stla(nst),freq(nfmax),sl(100),theta(100)
 integer Nfreq,Nsta,Nframes,Ifmin,Ifmax,cor

 integer narg
 narg=iargc()
 if (narg.ne.1) then
    write(*,*)'usage beamforming filelist param.dat'
    stop
 endif
!get filelist name
!-----------------
 call getarg (1,filein)
!get parameter from param.dat
!----------------------------
 open(9,file='param.dat')
 read(9,*) fmin,fmax,slmin,dsl,nsl,cor
 close(9)
 write(*,*) fmin,fmax,slmin,dsl,nsl,cor
! read file with data and botain spectrum 
!---------------------------------------------------
 beg=0.0
 open(10,file='filelist')
 do ista=1,nst ! loop over traces of NCF
    read(10,"(A300)",iostat=ferr)filein
    if (ferr.lt.0) then
       exit
    endif
!    write(*,*) filein
    call rsac1(filein,sig,npts,beg,delta,nsmax,nerr)
!    open(100,file='temp.dat')
!    write(100,"(F15.7)") (sig(i),i=1,npts)
!    close(100)
!    call getfhv('stla',stla(ista),nerr)
!    call getfhv('stlo',stlo(ista),nerr)
!    write(*,*) stlo(ista),stla(ista)
    if (nerr.ne.0) then
       write(*,*) 'Error reading in file: '
    endif
!    write(*,*) npts,beg,delta,nsmax,nerr
    !do FFT by FFTW for each NCF
    npow=ceiling((log(dble(npts))/log(2.0)));nsp=2**npow !determin power of FFT
    df=1/delta/nsp;
!    write(*,*) nsp,df
    do i=1,nsp/2+1
       freq(i)=(i-1)*df
    enddo
    do j=1,nsp 
       sigin(j)=sig(j)         
    enddo
!    open(11,file='temp.dat')
    call dfftw_plan_dft_r2c_1d(plan,nsp,sigin,sigout,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan,sigin,sigout)
    call dfftw_destroy_plan(plan)
    do j=1,nsp 
       fseis(j,ista)=sigout(j)
!      write(11,*) sigout(j)          
    enddo
!   close(11)

 enddo ! end loop over traces


!read x/y cooridinates
!----------------------------------------------------
 if(cor.eq.1) then
 open(11,file='xy.dat')
 do ista=1,nst
    read(11,*,iostat=ferr) stlo(ista),stla(ista)
    if (ferr.lt.0) then
       exit
    endif
 !   write(*,*) stlo(ista),stla(ista)
 enddo    
 endif
!calculate beamformer 
!----------------------------------------------------
 Nsta=ista-1
 Nfreq=nsp/2+1;
 Ifmin=int(fmin/df)+1
 Ifmax=int(fmax/df)+1

 do i=1,nsl
    sl(i)=slmin+(i-1)*dsl
!    write(*,*) sl(i)
 enddo
 do i=1,72
    theta(i)=(i-1)*5
!    write(*,*) theta(i)
 enddo
! do beamforming
 call fk(fseis,freq,stlo,stla,sl,theta,Nfreq,Nsta,Ifmin,Ifmax,nsl)
 

 return
 end

!========================================================
subroutine fk(fseis,freq,x,y,sl,theta,Nfreq,Nsta,Ifmin,Ifmax,nsl)

real ::  pi=3.14159265
complex fseis(32768,2000)
real freq(32768)
real, allocatable :: beam(:,:)
real x(100),y(100),sl(100),theta(100),beampower(72)
real projection(72,Nsta)
double complex rep(72,Nsta),temp1(72,Nsta),temp(72,72),repT(Nsta,72),cmat(Nsta,Nsta)
integer icc,iff,ith,itime,itc,Nfreq,Nsta,Nframes,Ntime,nsl
integer i,j,k,m,Ifmin,Ifmax
real isl,ifreq,omega
complex cj,tt
parameter(cj=(0.,1.))
allocate(beam(72,Ifmax-Ifmin+1))

write(*,*) Nfreq,Nsta,Ifmin,Ifmax

do i=1,72
      do j=1,Nsta
            projection(i,j)=x(j)*sin(theta(i)/180*pi)+y(j)*cos(theta(i)/180*pi);            
      enddo
enddo
!write(*,"(F15.7)") (projection(i,2),i=1,72)
open(100,file='beam.dat')

do icc=1,nsl
   isl=sl(icc)
   k=1
   do iff=Ifmin,Ifmax
      ifreq=freq(iff)  
      omega=ifreq*2*pi
      m=1;
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
         rep=exp((cj*omega*isl/1000)*projection) 
!         write(*,"(F15.7,F15.7)") (rep(1,i),i=1,Nsta)           
         repT=transpose(rep)
	 temp1=matmul(rep,cmat)
	 temp=matmul(temp1,conjg(repT))
         
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
