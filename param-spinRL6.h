!-------------------------------------------
!  num_proc=1-4  for energy minimization
!    use % mpirun if num_proc > 1
!-------------------------------------------
!
      integer*4   num_proc,np0,nbx,nhs,np6,mz6
      integer*4   meshx,meshy,meshz,P_max,MINTPOL
      real*8      av_start
!                          **
      parameter   (num_proc=6)
      parameter   (np0=12096,nbx=300,nhs=5000)   ! for  6x6x6
      parameter   (meshx=64,meshy=64,meshz=66)
      parameter   (np6=2016,mz6=11) ! np6= np0/4,mz6=16
      parameter   (P_max=3,MINTPOL=4*50048)
      parameter   (av_start=0.75d0)
!     parameter   (np0=8240,nbx=300,nhs=5000)    ! for  7x3x7
!     parameter   (np0=20330,nbx=300,nhs=5000)   ! for 11x3x11
!     parameter   (np0=113288,nbx=300,nhs=5000)  ! for 17x7x17
!
      character*4     praefixs*29,praefixi*29,praefixc*29, &
                      praefix*3,suffix*1
      common/filname/ praefixs,praefixi,praefixc,praefix,suffix
