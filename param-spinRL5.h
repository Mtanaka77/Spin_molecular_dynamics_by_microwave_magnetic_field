!  magnetite8.xyz, param-spinRL5.h 
!
      integer*4   num_proc,np0,nbx,nhs,meshx,meshy,meshz
      integer*4   P_max,MINTPOL
      real*8      av_start
!                          **
      parameter   (num_proc=5)
      parameter   (np0=7000,nbx=300,nhs=5000)    ! for 5x5x5
      parameter   (meshx=5,meshy=5,meshz=5)
      parameter   (P_max=3,MINTPOL=4*50048)
      parameter   (av_start=0.75d0)
!     parameter   (np0=8240,nbx=300,nhs=5000)    ! for  7x3x7
!     parameter   (np0=20330,nbx=300,nhs=5000)   ! for 11x3x11
!     parameter   (np0=113288,nbx=300,nhs=5000)  ! for 17x7x17
!
      character*4     praefixs*29,praefixi*29,praefixc*29,  &
                      praefix*3,suffix2*2,suffix1*2,suffix0*1
      common/filname/ praefixs,praefixi,praefixc,praefix,   &
                      suffix2,suffix1,suffix0
