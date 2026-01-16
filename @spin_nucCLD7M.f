c        1         2         3         4         5         6         71
c  Thread version
c    specify nthreads in param-spinRL.h
c
c  Large number of cells... 
c  1) Convergence is poor, so a small system is first organized. 
c   Then, a large system is equilibrated with the organized seed 
c   cell. 
c                                                   11/07/2010
c
c  2) Computation time becomes huge for the MC step, thus the sum 
c   is divided on processors, and they are mpi_allreduced.
c                                                   11/12/2010
c
c  3) Array arguments to /forces/ are passed by common statements,
c   as execution fails on mol205-208 (DL385G2).     11/27/2010
c
c
c  On SR16000
c  1. comment out - etime
c  2. in param.h for date_and_time, cdate*8, ctime*10 <- ctime*6 Linux
c  3. file directory: /msht/tanakam, /home/tanakam/MPI_spin (char*29)
c  4. Use n_smp=1 - thread version does not work!  Comment out
c     call fftw_threads...
c  5. use the script "f90sr3" with my /fftw3/lib
c
c  A) On SMP machines, call fftw_plan only on the major node,
c     Different nodes have different plan numbers.
c
******************************************************************
c*   Spin dynamics under the microwave H field for electrons,    *
c*   while nuclei feel forces from spins, neighboring nuclei     * 
c*   and restoring forces on their lattice points.               *          
c*   <- @spin_nuc.f.                                             *
c*** 6/27/2010 *********************************** 11/23/2010 ****
c*                                                               *
c*   Step 1: MC for equilibration at given T & B= 0.             *
c*   Step 2: MD (spin and nuclei dynamics)                       *
c*                                                               *
c*   If kstart= 0, MC is follwed by MD.                          *
c*      kstart > 0 using FT12 data                               *
c*---------------------------------------------------------------*
c*  >>NOTE for x()                                               *
c*                                                               *
c*       1  2 ... np1 ... np= np1+np2                            *
c*    1  Fe Fe Fe Fe  O   O                                      *
c*    2  Fe       Fe  O   O     rcut ... spins                   *
c*  ...  Fe       Fe  O   O     rcutC ... Coulomb                *
c*  np1  Fe Fe Fe Fe  O   O     rcutLJ... Lennard-Jones          *
c*  ...  O  O     O   O   O      - 3.7 Ang                       *
c*   np  O  O     O   O   O                                      *
c*                                                               *
c*  Distances                                                    *
c*   Fe(A)-O  1.891A, 3.495A                                     *
c*   Fe(B)-O  2.058A, 3.563A, 3.661A                             *
c*                                                               *
c*   Fe(A)-Fe(A)  3.636A                                         *
c*   Fe(B)-Fe(B)  2.969A                                         *
c*   Fe(A)-Fe(B)  3.481A                                         *
c*                                                               *
c*   O-O   2.85A, 2.97A, 3.088A                                  *  
c*                                                               *
c*    ... definition of i1(),i2() in spin_dynamics               *           
c*    ... multi write(06) must be avoided                        *
c*                                                               *
c* 1. Potential field derived to fit ab initio is used, in which *
c*   short-range force is A*exp(-Br) -C/r**6 +D/r**12            *  
c* 2. Fe(3+) and Fe(2+) are balanced (2+2) on each sub-lattice.  *
c*                                                               *
c************************************************ 02/11/2009 *****
c*                                                               *
c*  a) Give three axes of the crystal in /init/                  *
c*  b) If isolated crystal - do not fold in /spin_dynamics/      *
c*                                                               *
c*    Position of atoms: x,y,z (Ang)                             *
c*    Spin of atoms: spx,spy,spz (multiple of 1/2)               *
c*    Site index: site= 1,2 (prime number)                       *
c*     Positions of Fe(2+)/Fe(3+) at B-site are shuffled         *
c*                                                               *
c*    .............................                              *
c*               site=1   site=2                                 *
c*    .............................                              *
c*      spec=1   (A)3+    (B)3+                                  *
c*      spec=2            (B)2+                                  *
c*      spec=0             O 2-                                  *
c*    .............................                              *
c*                                                               *
c*   1. Energy minimization                                      *
c*     Apply the Metropolis criterion when Usys increases        *
c*     U= - Sum J_ij s_i*s_j +Sum_i g*mu_b B*s_i                 *
c*                                                               *
c*   2. Spin dynamics                                            *
c*     - should start from an energy-minimized state             *
c*                                                               *
c*    ds_i                 2*Jij       g*mue_b                   *
c*   ------ = s_i x (Sum_j ----- s_j - ------- Bex)              *
c*     dt                   hbar        hbar                     *
c*                                                               *
c*   For exchange interactions Jaa, Jbb and Jab, first and       *
c*  second neighbors ( r < 4A) are non-zero.                     *
c*                                                               *
c*  default units: 1 Ang, 10ps, 100 gauss                        *
c*                                                               *
c************************************************* 7/07/2007 *****
c
      include    'mpif.h'
      include    'param-spinRL.h'
c
      integer*4   COMM_LOCAL,size,rank,key,ierror,kstart,
     *            ir,ir0,igrp,ifbatch,wrt,n_MCsteps,nt_P3M
      real*4      cputime
      real*8      wtime,wtime1,wtime2
c
      common/parm5/  n_MCsteps,nt_P3M
      CHARACTER*8    LABEL,praefix*3,suffix*1
      COMMON/HEADR1/ LABEL,cdate
      common/print6/ ifbatch,wrt
      common/ranfff/ ir,ir0
      namelist/inp1/ kstart,praefix,suffix,ir0
c
c  kstart = 0 ... new run with MC (t=0), followed by MD. 
c         = 1 ... restart
c  
c  ---------------------------------------------
c  Switch  ft12 ---> ft10; copy to NFS volume /scratch6
c     kstart = 1
c  -------------------------------------------------
      label= 'sp3_nucl'
c
      call mpi_init (ierror)
      call mpi_comm_rank (MPI_COMM_WORLD,rank,ierror)
      call mpi_comm_size (MPI_COMM_WORLD,size,ierror)
c
c*******************************************
c*  A group for intra-task communication.  *
c*******************************************
c rank= 0... basic i/o
c
      if(size.eq.1) then
        write(06,*) 'stop: number of processors must be >= 2...'
        go to 900
      end if
c
      igrp = 1     ! Group #
      key= rank
      call mpi_comm_split (MPI_COMM_WORLD,igrp,key,COMM_LOCAL,ierror)
c
c -------------------------------------------------------------
c* Input by keyboard or input file < aaa
c
      if(rank.eq.0) then
c
        write(06,*) 'This run uses',size,' processors...'
        write(06,*) '         with',nthreads,' threads.'
c
        write(06,*) 'Type ifbatch (0,1) and kstart (0,1) ...'
        read(05,'(i1)') ifbatch   ! =0 to print by write(06)
        read(05,'(i1)') kstart 
c
        write(06,*) 'Give the file name by praefix and suffix'
        write(06,*) 'Type praefix (a3)...'
        read(05,'(a3)') praefix
c
        write(06,*) 'Type suffix (a1)...'
        read(05,'(a1)') suffix
c
        write(06,*) 'Type ir (i6)...'
        read(05,'(i5)') ir0
c
        write(06,*) 'kstart=',kstart
        write(06,*) 'praefix=',praefix
        write(06,*) ' Random seed: ir=',ir0
      end if
c
      call mpi_bcast ( kstart,1,mpi_integer,  0,MPI_COMM_WORLD,ierror)
      call mpi_bcast (praefix,3,mpi_character,0,MPI_COMM_WORLD,ierror)
      call mpi_bcast ( suffix,1,mpi_character,0,MPI_COMM_WORLD,ierror)
      call mpi_bcast ( ir0,1,mpi_integer,     0,MPI_COMM_WORLD,ierror)
      ir= ir0
c                                             + ++++++++++++++
      if(ifbatch.eq.0) then
        wrt= 06
      else
        wrt= 11
      end if
c
      praefixs = '/home/tanakam/MPI_spin/SAI'//praefix  ! NFS
      praefixi = '/msht/tanakam/heise7nc/sai'//praefix   ! NFS
      praefixc = '/msht/tanakam/heise7nc/sai'//praefix
c -------------------------------------------------------------
c
      nt_P3M = 5
c
      if(rank.eq.0) then
        if(ifbatch.eq.1) OPEN (unit=wrt,file=praefixc//'.06'//suffix,
     *                             status='unknown',form='formatted')
c
        call date_and_time (cdate,ctime)
        write(wrt,*) ' This run uses ',size,' processors...'
        write(wrt,*) '         with',nthreads,' threads.'
        write(wrt,*) '         '
        write(wrt,*) '   P3M is called in every',nt_P3M,' steps'
        write(wrt,*) '         '
        write(wrt,*) '   kstart =',kstart
        write(wrt,*) '   ifbatch =',ifbatch
        write(wrt,*) '         '
        write(wrt,*) '   today = ',cdate
        write(wrt,*) '   time  = ',ctime
c
        if(ifbatch.eq.1) close (wrt)
c
        OPEN (unit=77,file=praefixc//'.77'//suffix,
     *                             status='unknown',form='formatted')
c
        nframe= 4
        call gopen (nframe)
        close (77)
      end if
c
      call clocks (cputime,wtime1)
c ------------------------------------------------------------------
      call spin_dynamics (rank,size,igrp,COMM_LOCAL,kstart,suffix)
c ------------------------------------------------------------------
      call clocks (cputime,wtime2)
c
      if(rank.eq.0) then
        if(ifbatch.eq.1) OPEN (unit=wrt,file=praefixc//'.06'//suffix,
     *            status='old',position='append',form='formatted')
c
        wtime = wtime2 -wtime1
        write(wrt,693) cputime,wtime
  693   format(/,' cpu, wtime (sec)=',2f10.3)
c
        if(ifbatch.eq.1) close (wrt)
      end if
c*
  900 continue
      call mpi_comm_free (COMM_LOCAL,ierror)
      call mpi_finalize  (ierror)
c
      stop
      end
c
c
c------------------------------------------------------------------------
      subroutine spin_dynamics (rank,size,igrp,COMM_LOCAL,kstart,suffix)
c------------------------------------------------------------------------
      implicit none
c
      include    'mpif.h'
      include    'param-spinRL.h'
c
      integer*4  COMM_LOCAL,rank,size,igrp,kstart,np1,np2,
     *           i1(0:num_proc),i2(0:num_proc),i3(0:num_proc),
     *           i4(0:num_proc),
     *           cnt_recv(0:num_proc),disp_recv(0:num_proc),
     *           cnt_recvC(0:num_proc),disp_recvC(0:num_proc),
     *           cnt_send,ierror,MPIerror,i_mc,ifcmp,kk
c
      real*8     x,y,z,spx,spy,spz,ch,fx,fy,fz,fxC,fyC,fzC,
     *           Jint,r_ij,rintC,fc1,fc2,fcLJ
      common/partcl1/ 
     *           x(np0),y(np0),z(np0),spx(np0),spy(np0),spz(np0),
     *           ch(np0)
      common/partcl2/ 
     *           fx(np0),fy(np0),fz(np0),fxC(np0),fyC(np0),
     *           fzC(np0),Jint(nbx,np0),r_ij(nbx,np0),rintC(nbx,np0)
c
      real*8     vx(np0),vy(np0),vz(np0),mass(np0),ag(np0),ep(np0),
     *           fkx(np0),fky(np0),fkz(np0),sp2(np0),spx0(np0),
     *           spy0(np0),spz0(np0),spx1(np0),spy1(np0),spz1(np0),
     *           wsp1(3,np0),wsp2(3,np0),Jint_1,
     *           U,U1,U2,spx00(np0),spy00(np0),spz00(np0),
     *           U_1,U_2,aspx(np0),aspy(np0),aspz(np0)
c
      integer*4  nintS,lintS,nintC,lintC,if_LJ,spec,site
      common/partcl3/ 
     *           nintS(np0),lintS(nbx,np0),nintC(np0),lintC(nbx,np0),
     *           if_LJ(np0),spec(np0),site(np0)
c
      integer*4  nlist(np0),lmax,k0,iac1,iac2,irej,modes,kwrite,
     *           n_of_J,np10,np100,ia,ja,ka,i1x,i2x,i1y,i2y,i1z,i2z,
     *           n_MCsteps,nt_P3M
c
      real*8     t,dt,dts,dt0,dth,bex,bey,bez,spin2,spin3,
     *           Jaa,Jbb,Jab,J00,B00,Bapx,Bapy,Bapz,bap,
     *           tau_b,tau_R,tau_diss,t_adv,fw,fw00,
     *           Temp,Tcurie,g,mue_b,hbar,kJoule,kcal,mol,kT,eV,
     *           omg_b,rnd,
     *           f1,g1,b1,prb,qsx,qsy,qsz,rsx,rsy,rsz,hh1,hh2,
     *           qq,rqq,dthg,dthqq,rdiv,spmin,
     *           alp,xx,yy,zz,rr,rr0,rcut,toler,deps,Usys8(nhs),
     *           ss,smax,sav,t_unit,e_unit,a_unit,m_unit,m_fe,m_o,
     *           pi,pi2,th,ph,tmax,tmax0,cptot,unif1(2),unif2(2),
     *           xleng0,yleng0,zleng0,sss,rlist(np0),
     *           buffer1(4),buffer2(4),del_en,wtime,
     *           fc3,J_ki,wfdt,vth0,vth_o,vth_f,vth,
     *           svx,svy,svz,sqrt2,vmax1,dgaus2,vsq1,vsq2,
     *           vx1,vy1,vz1,mas,E_sp,E_C_r,E_LJ,E_Coulomb_P3M
c
      real*8     x0(np0),y0(np0),z0(np0),vx0(np0),vy0(np0),vz0(np0)
      real*8     x_0(np0),y_0(np0),z_0(np0)
      real*4     dtwr,dtwr2,cputime
c
      real*8        rad_fe,rad_o,elj_fe,elj_o,rcutLJ,rcutC
      common/atoms/ rad_fe,rad_o,elj_fe,elj_o,rcutLJ,rcutC
c
      integer*4     it,is,iw,iwa,iwb,istop,ifbatch,wrt,
     *              iwrt1,iwrt2,iwrta,iwrtb,iter,itermax,
     *              i00,i,j,k,l,nsite,notconv,nframe,mx,my,mz,
     *              ir,ir0,lsite(num_proc,np0),ic,nstep_mc,nstep_mcT,
     *              ifdt,kmax,if_vres(np0),ifv
      common/parm1/ it,is
      common/parm2/ dtwr,dtwr2
      common/parm3/ t,rcut,pi,dt,tmax,cptot
      common/parm4/ mx,my,mz
      common/parm5/ n_MCsteps,nt_P3M
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,
     *              tau_diss,Temp,Tcurie
      common/imemo/ iwa,iwb
      common/ranfff/ ir,ir0
      common/itera/  toler,itermax
      common/print6/ ifbatch,wrt
c
      real*4        spinx,spinz,spin7,bextx,bextz,magx,magy,magz,
     *              Usys,conv,aitr,psdt,Tfix,Uss,Usb,Tsx,Tsy,Tsz,sum_mb,
     *              U_fe,U_o,ds_fe,ds_o,fdt4,vdt4,idt4,timeh,dtrhs
      real*4        sx1(3),sy1(3),sz1(3),sx2(3),sy2(3),sz2(3),sn1(3),
     *              csx,csy,csz,axis(100),freq(100),t4,bex4,bez4,
     *              spx4(np0),spy4(np0),spz4(np0),
     *              ch4(np0),x4(np0),y4(np0),z4(np0),
     *              vx4(np0),vy4(np0),vz4(np0)
      real*8        wx,wy,wz,wn,wsx,wsy,wsz,wp,
     *              wx1,wy1,wz1,wn1,Uav,wt1,wx7,wy7,wz7,wn7,Uav7,wt7,
     *              ssx,ssy,ssz,tht,phi,psi,bbw,Tsz0,aTsz,av_tsz(nhs),
     *              Ub,Um,Ub1,Um1,UU1(3),UU2(3),ranff,fdt8,vdt8,ds1,ds2
      common/ehist/ spinx(nhs),spinz(nhs),spin7(nhs),
     *              bextx(nhs),bextz(nhs),magx(nhs),magy(nhs),magz(nhs),
     *              Usys(nhs),conv(nhs),aitr(nhs),psdt(nhs),Tfix(nhs),
     *              Uss(nhs),Usb(nhs),Tsx(nhs,3),Tsy(nhs,3),Tsz(nhs,3),
     *              sum_mb(nhs),U_fe(nhs),U_o(nhs),ds_fe(nhs),ds_o(nhs),
     *              fdt4(nhs),vdt4(nhs),idt4(nhs),timeh(nhs)
      logical       MC_first
      data          MC_first/.true./
c
      real*8         alpha,xleng,yleng,zleng
      integer*4      PP
      common/ewald1/ alpha,fc2
      common/ewald2/ PP
      common/ewald3/ xleng,yleng,zleng
c
      character*8     praefix*3,suffix*1,LABEL,
     *                plot_ch*8,char*2
      real*4          time
      COMMON/HEADR1/  LABEL,cdate
      COMMON/HEADR2/  time
      namelist/inp1/  kstart,praefix,suffix,ir0
c
      istop = 0
c
      if(rank.eq.0) then 
        if(ifbatch.eq.1) OPEN (unit=wrt,file=praefixc//'.06'//suffix,
     *               status='old',position='append',form='formatted')
c
        if(kstart.eq.0) then
          write(wrt,21) 
   21     format(/,' ## MC step with B_mw=0, followed by MD.....',/)
        else
          write(wrt,23) 
   23     format(/,' ## MD step by restart.....',/)
        end if
      end if
c      
c--------------------------
c*  Cold start: MC and MD
c--------------------------
      if(kstart.eq.0) then
        t= 0.d0
c
        it= 0
        is= 0   ! write out history
c
        iwa= -1
        iwb= -1
c
c Parameters are read
c
        call READ_CONF (xleng0,yleng0,zleng0,rank,suffix)
c
c xleng is expanded and redeined in "init"
c
        call init (x,y,z,ch,spx,spy,spz,sp2,spec,site,
     *             xleng0,yleng0,zleng0,np1,np2,np10,
     *             i1x,i2x,i1y,i2y,i1z,i2z,rank,size,suffix)
c
        do i= 1,np1+np2
        x_0(i)= x(i)
        y_0(i)= y(i)
        z_0(i)= z(i)
        end do
c
        if(np1+np2.gt.np0) then
          write(wrt,*) ' # Stop: np > np0...'
          return
        end if
c
c
c--------------------------
c*  Restart from FT10
c--------------------------
      else
        if(rank.eq.0) write(wrt,*) 'Read: ',praefixi//'.10'//suffix
c
        OPEN (unit=10,file=praefixi//'.10'//suffix,
     *                             status='old',form='unformatted')
c                  +++++ +++++ +++++
        read(10) t,xleng,yleng,zleng,rcut,rcutC,Temp,Tcurie,
     *           tmax,dt,cptot
        read(10) x,y,z,vx,vy,vz,ch,mass,ag,ep
        read(10) x_0,y_0,z_0,rintC
        read(10) spx,spy,spz,sp2,spx00,spy00,spz00,r_ij
        read(10) aspx,aspy,aspz
        read(10) Jint,U,spin2,spin3
        read(10) Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,toler,itermax
        read(10) dtwr,dtwr2,fw00,aTsz,Tsz0,av_tsz
        read(10) it,is,iwa,iwb,ir,ir0
        read(10) np1,np2,spec,site,nintS,lintS,nintC,lintC,if_LJ
        read(10) i1,i2,i3,i4,disp_recv,cnt_recv,disp_recvC,cnt_recvC
        read(10) spinx,spinz,spin7,bextx,bextz,magx,magy,magz,
     *           Usys,conv,aitr,psdt,Tfix,Uss,Usb,Tsx,Tsy,Tsz,sum_mb,
     *           U_fe,U_o,fdt4,vdt4,idt4,timeh
        close(10)
c
c---------------------------
c* Overwrite above values
c---------------------------
        call READ_CONF (xleng0,yleng0,zleng0,rank,suffix)
      end if
c
c
      if(rank.eq.0) then
        if(kstart.eq.1) write(wrt,*) ' Spin dynamics ................'
c
        write(wrt,*) '>>Input parameters...'
        write(wrt,931) mx,my,mz
        write(wrt,932) kstart
        write(wrt,933) Bapx,Bapy,Bapz,Temp,
     *              10*Jaa,10*Jbb,10*Jab,rcut,rcutC
        write(wrt,934) dt,toler
  931   format(' Domains: mx, my, mz=',3i3)
  932   format(' kstart=',i3,/)
  933   format(' Applied B field(/100gauss) = ',3f8.1,/,
     *         ' Temp(K) = ',0pf7.1,/,
     *         ' Jaa(meV) = ',1pe11.3,/,
     *         ' Jbb(meV) = ',e11.3,/,
     *         ' Jab(meV) = ',e11.3,/,
     *         ' rcut(Ang)= ',e11.3,/,
     *         ' rcut_Coul(Ang)= ',e11.3,/)
  934   format(' dt  = ',  1pe11.3,/,
     *         ' tolerance = ',1pe11.3,/)
      end if
c
c--------------------------
c* Define constants
c--------------------------
      pi = 4*datan(1.d0)
      pi2= 2*pi
c
      t_unit= 1.d-12         ! ps
      a_unit= 1.d-8          ! Ang
      e_unit= 4.80325d-10    ! e
      m_unit= 1.67261d-24    ! m_H
c
      hbar = 6.62620d-27/(2*pi)  ! erg*sec
      mue_b= 9.27410d-21         ! e*hbar/2mc
      kJoule  = 1.d10            ! erg
      kcal= 4.1868d0 *kJoule     ! 4.18 J/cal
      mol = 6.0220d23
c
      g  = 2.0023d0
      kT = 1.38062d-16*Temp
      eV = e_unit/300.d0
c
      m_fe= 56
      m_o = 16
c
      vth0= (t_unit/a_unit) * sqrt(2*kT/m_unit)  ! for Hydrogen atom
      vth_o= vth0/sqrt(m_o)
      vth_f= vth0/sqrt(m_fe)
c
c* Exchange: J= kT_c/0.3z   cf. Kittel Eq.(15-7)
c
c     znei= 6
c     asp = sp2(1)
c     J00= 3.0d0* (3/(2*znei*asp*(asp+1))) *1.38062d-16*Tcurie
      J00 = 1.0d-2 ! 10 meV
        if(rank.eq.0) write(wrt,*) ' J00(eV)=',J00
c
      J00 = J00 * eV
      B00 = 100.d0  ! unit: 100 gauss
c
      f1 = J00 *t_unit/hbar
      g1 = g*mue_b*B00 *t_unit/hbar
c
      bap= sqrt(Bapx**2 +Bapy**2 +Bapz**2)
      b1= g*mue_b*B00 *bap  ! mue*B energy
c
      tau_R= 0.
      if(bap.ne.0) tau_R= 1.d12* pi2/((g*mue_b*(B00*bap))/hbar)  ! ps
      omg_b= pi2/tau_b
c
      if(rank.eq.0) then
        write(wrt,101) Tcurie,J00
        write(wrt,102) (J00*Jab)*spin2**2,b1*spin2,kT,omg_b
  101   format(' Tcurie=',f7.1,/,
     *         ' J00       = ',1pe12.5)
  102   format(' Jab*s*s   = ',1pe12.5,/,
     *         ' g*mub*s*B = ',  e12.5,/,
     *         ' kT(erg)   = ',  e12.5,/,
     *         ' omg_b(/ps)= ',  e12.5,/)
c
        dtrhs= f1*Jab*sp2(1)**2*dt
        if(abs(dtrhs).gt.0.4) then
          write(wrt,*) ' ## dt*rhs > 0.4 ..... stop'
          return
        end if
c
        write(wrt,103) dtrhs,tau_R,tau_b,tau_diss
  103   format(' Increment of RHS for dt=',f10.4,
     *         /,'   tau_res = ',f10.2,
     *         /,'   tau_mw  = ',f10.2,
     *         /,'   tau_diss= ',f10.2,/)
c
        write(wrt,104) f1*Jab*sp2(1)**2,
     *                g1*bap*sp2(1)
  104   format('## RHS of spin equation (per t_unit/hbar)...',/,
     *         '   J*si*sj = ',f10.6,/,
     *         ' g*mue*B*s = ',f10.6,/)
c
        write(wrt,105) vth_o
  105   format('## Thermal velocity of Oxygen (a_unit/t_unit)...',/,
     *         '   v_th(O)=',1pe12.3,/)
      end if
c
c ----------------------------
c*  Force constants
c ----------------------------
c ############################
c  Unit: Coulomb force ->(e*tau/a)^2/m
c
      fc1 = J00*(t_unit/a_unit)**2/m_unit          ! spin force
      fc2 = (e_unit*t_unit)**2/(m_unit*a_unit**3)  ! Coulomb force
      fcLJ = 48.d0*t_unit**2*kT/(m_unit*a_unit**2)
c
      fc3= m_unit*(a_unit/t_unit)**2    ! for kinetic energy
c
c* Order of J(r) as (1/r)**n
c
      n_of_J= 1
c
      if(rank.eq.0) then
        write(wrt,107) fc1,fc2,fc3,n_of_J
  107   format('## Force constants on nuclei ###',/,
     *         '  fc_spin=',1pe12.4,/
     *         '  fc_Coul=',e12.4,/
     *         '  fc_kin =',e12.4,/
     *         '  n_of_J(r)=',i3,/)
      end if
c ############################
c
c  >> P3M_init ...
c* Define parameters of /ewald1-3/
c  -------------------------------
      pi = 4.d0*DATAN(1.d0)
cc    alpha  = 2*pi/zleng   ! 0.26d0
      alpha  = 2*pi/rcutC   ! near if r < rcutC ?
      PP     = P_max 
c  -------------------------------
c
c   fc2 -> prefactor= (t_unit*e_unit)**2/(w_unit*a_unit**3)
c     pref2    = 48.d0*t_unit**2*kbT/(m_unit*a_unit**2)
c
      call  interpol_charge_assign_function (rank,wrt)
      call  calculate_meshift (rank,wrt)
      call  calculate_differential_operator (rank,wrt)
      call  calculate_influence_function (rank,wrt)
c
      if(rank.eq.0) then
        write(wrt,601) rcutC,rcutLJ,mx,my,mz 
  601   format('## P3M parameters ###',/,
     *         '  rcutC  =',f7.2,/
     *         '  rcutLJ =',f7.2,/
     *         '  [mx, my, mz] domains =',3i3,/)
c
        write(wrt,*) ' P3M successfully initialized !'
      end if
c
c
c--------------------------
c* Parallelization index
c--------------------------
c                  **
      do k= 0, size-1
      i1(k)=     k*(np1/size) +1
      i2(k)= (k+1)*(np1/size)
      if(k.eq.size-1) i2(k)= np1
c
      disp_recv(k)= 3*(i1(k) -1)
      cnt_recv(k) = 3*(i2(k) -i1(k) +1)
c
      if(kstart.eq.0 .and. rank.eq.0) then
        if(k.eq.0) then
          write(wrt,*) '  '
          write(wrt,*) 'k, i1,i2, disp_recv, cnt_recv...'
        end if
        write(wrt,29) k,i1(k),i2(k),disp_recv(k),cnt_recv(k)
   29   format('   k=',i3,4i8)
      end if
      end do
c
c                  **
      do k= 0, size-1
      i3(k)=     k*((np1+np2)/size) +1
      i4(k)= (k+1)*((np1+np2)/size)
      if(k.eq.size-1) i4(k)= np1 +np2
c
      disp_recvC(k)= 3*(i3(k) -1)
      cnt_recvC(k) = 3*(i4(k) -i3(k) +1)
c
      if(kstart.eq.0 .and. rank.eq.0) then
        if(k.eq.0) then
          write(wrt,*) '  '
          write(wrt,*) 'k, i3,i4, disp_recvC, cnt_recvC...'
        end if
        write(wrt,29) k,i3(k),i4(k),disp_recvC(k),cnt_recvC(k)
      end if
      end do
c
c
c------------------------
c* Interaction table (1)
c------------------------
c Pick up only nearest neighbors for A-A, B-B, A-B
c depending on my spin
c   site(i)= 1   ! A site (8)
c   site(i)= 2   ! B site (16)
c
      do 30 i= 1,np1
      l= 0
c
      do 32 j= 1,np1
      if(j.eq.i) go to 32   ! avoid itself
c
      xx= x_0(i) -x_0(j)
      yy= y_0(i) -y_0(j)
      zz= z_0(i) -z_0(j)
c
c* Do not fold - if isolated crystal
c  #################################
c
      xx= xx -DNINT(xx/xleng)*xleng
      yy= yy -DNINT(yy/yleng)*yleng
      zz= zz -DNINT(zz/zleng)*zleng
      rr= sqrt(xx**2 +yy**2 +zz**2)
c
      if(rr.gt.rcut) go to 32
      l= l +1
      rlist(l)= rr
      nlist(l)= j
   32 continue
c
      lmax= l
c
      l= 0
   40 l= l +1
      if(l.gt.lmax) go to 30
c
      rr0= 1000
      if(l.gt.nbx) then
        write(wrt,*) 'l > nbx ... stop at 40'
        return
      end if
c
c Sort the indices in ascending order
      do 50 k= 1,lmax
      if(rlist(k).le.rr0) then
        rr0= rlist(k)
        j  = nlist(k)
        k0 = k
      end if
   50 continue
c
      nintS(i)= l
      lintS(l,i)= j
      r_ij(l,i)= rr0
c
c  Exchange interaction is small for distant pairs: rr > 3 Ang
c
      nsite= site(i)*site(j)
      lsite(l,i)= nsite
      if(nsite.eq.1) Jint(l,i)= Jaa
      if(nsite.eq.2) Jint(l,i)= Jab
      if(nsite.eq.4) Jint(l,i)= Jbb
c
      rlist(k0)= 1000       ! exclude in the next search
      go to 40
   30 continue
c
c------------------------
c* Interaction table (2)
c------------------------
c  Coulomb force
c
c    Fe(B)  2.058 x 6 - O  *
c           2.969 x 6 - Fe *
c           3.481 x 6 - Fe
c
c    Fe(A)  1.891 x 4  - O *
c           3.481 x 12 - Fe
c           3.495 x 12 - O
c
c      O    1.891 x 1 - Fe *
c           2.058 x 3 - Fe *
c           2.850 x 3 - O  *
c           2.970 x 6 - O  *
c           3.088 x 3 - O  *
c           3.495 x 3 - Fe    
c
c  Give the LJ-pair limit in SAI100 data file
c    r_ij(t=0) < rcutLJ = rcutC
c
c* Divided upon processors for a large system... (11/12/2010)
c
      do 60 i= i3(rank),i4(rank)
      l= 0
c
      do 62 j= 1,np1+np2
      if(j.eq.i) go to 62   ! avoid itself
c
      xx= x_0(i) -x_0(j)
      yy= y_0(i) -y_0(j)
      zz= z_0(i) -z_0(j)
c
      xx= xx -DNINT(xx/xleng)*xleng
      yy= yy -DNINT(yy/yleng)*yleng
      zz= zz -DNINT(zz/zleng)*zleng
      rr= sqrt(xx**2 +yy**2 +zz**2)
c
      if(rr.gt.rcutC) go to 62
      l= l +1
      rlist(l)= rr
      nlist(l)= j
   62 continue
c
      lmax= l
c
      l= 0
   70 l= l +1
      if(l.gt.lmax) go to 60
c
      if(l.gt.nbx) then
        write(wrt,*) 'l > nbx ... stop at 70, nbx=',nbx
        return
      end if
      rr0= 1000
c
c Sort the indices in ascending order
      do 80 k= 1,lmax
      if(rlist(k).le.rr0) then
        rr0= rlist(k)
        j  = nlist(k)
        k0 = k
      end if
   80 continue
c
      nintC(i)= l
      lintC(l,i)= j
      rintC(l,i)= rr0
c
c  Use the same list, but limit the LJ pair 
c     for r_ij(0) < 3.7 Ang (almost 2nd neighbors) -> 7.0 Ang for 5x5x5 cells
c                                    9/18/2010                10/27/2010
      if_LJ(i)= 0
      if(rr0.lt.rcutLJ) if_LJ(i)= 1
c
c
      rlist(k0)= 1000       ! exclude in the next search
      go to 70
   60 continue
c
c
c%%%% only at t=0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%
      if(kstart.eq.0) then
c
c ------------------------------------------------
c*  Mass, radius, LJ, and initial velocity
c ------------------------------------------------
c
      do i= 1,np1+np2
      if(i.le.np1) then
        mass(i)= m_fe
        ag(i)= rad_fe
        ep(i)= elj_fe*(kcal/mol)/kT
      else
        mass(i)= m_o
        ag(i)= rad_o
        ep(i)= elj_o*(kcal/mol)/kT
      end if
      end do
c
c***************
c*  Velocity.  *
c***************
c
      sqrt2= dsqrt(2.d0)
      call ggauss
c
c* For irons
c
      svx= 0.
      svy= 0.
      svz= 0.
c
      do i= 1,np1
      vmax1= sqrt2*vth0/sqrt(mass(i))
      vx(i)= dgaus2(vmax1)
      vy(i)= dgaus2(vmax1)
      vz(i)= dgaus2(vmax1)
c
      svx= svx +vx(i)
      svy= svy +vy(i)
      svz= svz +vz(i)
      end do
c
      do i= 1,np1
      vx(i)= vx(i) -svx/np1
      vy(i)= vy(i) -svy/np1
      vz(i)= vz(i) -svz/np1
      end do
c
c* For oxygens
c
      svx= 0.
      svy= 0.
      svz= 0.
c
      do i= np1+1,np1+np2
      vmax1= sqrt2*vth0/sqrt(mass(i))
      vx(i)= dgaus2(vmax1)
      vy(i)= dgaus2(vmax1)
      vz(i)= dgaus2(vmax1)
c
      svx= svx +vx(i)
      svy= svy +vy(i)
      svz= svz +vz(i)
      end do
c
      do i= np1+1,np1+np2
      vx(i)= vx(i) -svx/np2
      vy(i)= vy(i) -svy/np2
      vz(i)= vz(i) -svz/np2
      end do
      end if  
c%
c%%%% only at t=0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-------------------------------------------------------
      if(rank.eq.0) then
        if(ifbatch.eq.1) close (wrt)
c
        write(21,*) '### Spin tables...'
        do 901 i= 1,np1
        write(21,936) i,sp2(i),site(i),nintS(i)
        write(21,937) (lintS(k,i),k=1,nintS(i))
        write(21,938) (lsite(k,i),k=1,nintS(i))
        write(21,939) (r_ij(k,i),k=1,nintS(i))
  936   format('i=',i5,/,'  spin,site,nint=',f5.1,2i3)
  937   format(' lintS= ',15i5)
  938   format(' lsite=',15i5)
  939   format(' r_ij=',10f6.2)
  901   continue
c
        write(21,*) '  '
        write(21,*) '### Coulomb force pairs...'
        do 910 i= 1,np1+np2
        if(i.eq.1) write(21,*) '>> For Iron ...'
        if(i.eq.np1+1) write(21,*) '>> For Oxygen ...'
c
        if(i.le.np1) then
          if(site(i).eq.1) char= ' A'
          if(site(i).eq.2) char= ' B'
        else
          char= '  '
        end if
c
        write(21,940) i,nintC(i),char,site(i),spec(i)
        write(21,941) (k,lintC(k,i),rintC(k,i),k=1,nintC(i))
  940   format('i=',i5,' nintC=',i3,' site(A,B)=',a2,' site, spec=',2i2)
  941   format(' k,j,r_ij= ',i3,i5,f10.3)
  910   continue
c
        OPEN (unit=13,file=praefixc//'.13'//suffix,
     *                              status='unknown',form='unformatted')
c
        write(13) xleng,yleng,zleng,rcut,rcutC,dt
        write(13) Jaa,Jbb,Jab,Bapx,Bapy,Bapz,Jint
        write(13) np1,np2,spec,site,nintS,lintS
        write(13) nintC,lintC,if_LJ
        close(13)
      end if
c-------------------------------------------------------
c-------------------------------------
c* Time or minimization loop begins
c-------------------------------------
c
      if(kstart.eq.1) del_en= 0
c
 1000 continue
      dth= 0.5d0*dt
c
      it= it + 1
      t = t + dt
c
c* All nodes must share the same info
c  Time synchronization is achieved here...
c
      call clocks (cputime,wtime)
      if(t.gt.tmax) istop= 1
      if(wtime/3600.d0.gt.cptot) istop= 1
c
      call mpi_bcast (istop,1,mpi_integer,0,MPI_COMM_WORLD,ierror)
      if(istop.eq.1) go to 7000
c
c     if(kstart.ne.0) then
c       dtwr =  tau_b/80.d0
c       dtwr2 = tau_b/4.d0
c     end if
c
      iwrt1= iwrta(t,dtwr)
      iwrt2= iwrtb(t,dtwr2)
      iter = 0
c
      if(it.eq.1) then
        wx1= 0
        wy1= 0
        wz1= 0
        wn1= 0
        Uav= 0
        wt1= 0
      end if
c
c
c----------------------------------------------------------------
c* MC step is executed by all nodes if kstart= 0
c   Arrays spx00...spz00() are used in MC - do not use spx-spz !
c----------------------------------------------------------------
c
 2100 continue
      if(kstart.eq.0) then
c
        if(MC_first) then   ! Organize a seed one cell
          np100= np10       ! one cell
          nstep_mc = 50001
          kwrite   = 10000
        else                ! Initialize a large system with the organized seed cell
          np100= np1        ! all
          nstep_mc = n_MCsteps  ! 1000001
          kwrite   =   50000    ! 50000
        end if
c
        if(rank.eq.0) then
          if(ifbatch.eq.1) OPEN (unit=wrt,file=praefixc//'.06'//suffix,
     *                  status='old',position='append',form='formatted')
          if(MC_first) then
            write(wrt,*) 'Step 1: MC for a small (one cell) system...'
          else
            write(wrt,*) '   '
            write(wrt,*) 'Step 2: Equilibration of a large system...'
            write(wrt,*) '     maximum iteration count=',n_MCsteps
          end if
          if(ifbatch.eq.1) close (wrt)
        end if
c
        do i= 1,np1
        spx00(i)= spx(i)
        spy00(i)= spy(i)
        spz00(i)= spz(i)
        end do
c
        go to 2000
      else
c
        if(it.eq.1) then
          call magnetiz (spx,spy,spz,g,wx7,wy7,wz7,wn7,U1,Uav7,
     *                   wt7,np1,-1)
c
c* Sign of magnetic field is defined here
c
          if(wz7.gt.0) then
            fw00= 1.d0    ! choose so that -M*B < 0
          else
            fw00= -1.d0
          end if
        end if
c
        go to 3000
      end if
c
c--------------------------------------
c* Energy minimization by MC
c--------------------------------------
c  Flip a spin
c
 2000 i_mc= 0
c
      bex= 0
      bey= 0
      bez= 0
c
      U= 1.d10
      iac1= 0
      iac2= 0
      irej= 0
c
 2300 i_mc= i_mc +1
      j= mod(int(np100*ranff(0.)+0.00001),np100) +1
c                  +++                      +++
      th= pi*ranff(0.)
      ph= pi2*ranff(0.)
c
      spx0(j)= spx00(j)   ! save
      spy0(j)= spy00(j)
      spz0(j)= spz00(j)
c
      spx00(j)= sp2(j)*dsin(th)*dcos(ph)
      spy00(j)= sp2(j)*dsin(th)*dsin(ph)
      spz00(j)= sp2(j)*dcos(th)
c
      U1= 0.d0
c
c* Calculating pairs both ways - 0.5*
c
      if(MC_first) then
c
        do 200 i= 1,np100
        U1= U1 + b1 * (bex*spx00(i) +bey*spy00(i) +bez*spz00(i))
c
        do 210 l= 1,np10
        if(l.eq.i) go to 210
c
        xx= x(i)-x(l)
        yy= y(i)-y(l)
        zz= z(i)-z(l)
c                             +       + small system
        xx= xx -DNINT(xx/xleng0)*xleng0
        yy= yy -DNINT(yy/yleng0)*yleng0
        zz= zz -DNINT(zz/zleng0)*zleng0
        rr= sqrt(xx**2 +yy**2 +zz**2)
c
        if(rr.gt.rcut) go to 210
        nsite= site(i)*site(l)
c
        if(nsite.eq.1) Jint_1= Jaa
        if(nsite.eq.2) Jint_1= Jab
        if(nsite.eq.4) Jint_1= Jbb
        U1= U1 - 0.5d0*J00* Jint_1*(spx00(i)*spx00(l) 
     *                        +spy00(i)*spy00(l) +spz00(i)*spz00(l))
  210   continue
  200   continue
c
      else
c
c* Divided upon processors for a large system... (11/12/2010)
c
        do 220 i= i1(rank),i2(rank)
        U1= U1 + b1 * (bex*spx00(i) +bey*spy00(i) +bez*spz00(i))
c
        do 220 k= 1,nintS(i)
        l= lintS(k,i)
        U1= U1 - 0.5d0*J00* Jint(k,i)*(spx00(i)*spx00(l) 
     *                        +spy00(i)*spy00(l) +spz00(i)*spz00(l))
  220   continue
c
c
        U_1= U1
        call mpi_allreduce (U_1,U_2,1,mpi_real8,mpi_sum,COMM_LOCAL,
     *                      MPIerror)
        U1= U_2
      end if
c
c
c* Apply the Metropolis criterion when U increases
c
      if(U1.lt.U) then            ! Accept the flip
        U= U1                     ! Keep lowest energy
        iac1= iac1 +1
      else
c
        prb= exp(-(U1-U)/kT)       ! (U1-U) for one excitation
        if(prb.gt.ranff(0.)) then
          iac2= iac2 +1
          U= U1                   ! <-- Keep this energy
        else
          irej= irej +1           ! Reject the trial flip and
          spx00(j)= spx0(j)         ! restore the spin
          spy00(j)= spy0(j)
          spz00(j)= spz0(j)
        end if
      end if
c
      if(rank.eq.0) then
        call magnetiz (spx00,spy00,spz00,g,wx7,wy7,wz7,wn7,U1,Uav7,
     *                 wt7,np100,1)
c                            +++
        if(mod(i_mc,kwrite).eq.1) then
          if(ifbatch.eq.1) OPEN (unit=wrt,file=praefixc//'.06'//suffix,
     *               status='old',position='append',form='formatted')
c
          call magnetiz (spx00,spy00,spz00,g,wx7,wy7,wz7,wn7,U1,Uav7,
     *                   wt7,np100,7)
c                              +++
          bbw= B00*bez
          write(wrt,201) i_mc,iac1,iac2,irej,U/np100,U1/np100,
     *                  bbw,wx7,wy7,wz7     !    +++      +++
  201     format(' i_mc,iac1,iac2,irej=',i7,3i7,' <U>,<U1>=',
     *           1p2e11.3,' Bz=',e9.2,' <m>=',0p3f6.2)
c
          call magnetiz (spx00,spy00,spz00,g,wx7,wy7,wz7,wn7,U1,Uav7,
     *                   wt7,np100,0)
c                              +++
          if(ifbatch.eq.1) close (wrt)
        end if
      end if
c
      if(i_mc.lt.nstep_mc) go to 2300
c     -------------------------------
c
      if(MC_first) then
        MC_first= .false.
c
        i= 0
        do 230 ia= i1x,i2x
        do 230 ja= i1y,i2y
        do 230 ka= i1z,i2z
c
c To enable convergence to equilibrium...  11/13/2010
c
        if(mod(ia,3).eq.0 .and.
     *     mod(ja,3).eq.0 .and.
     *     mod(ka,3).eq.0) then
c
c Seed cell
          do 240 l= 1,np10
          i= i +1
          spx(i)= spx00(l)  ! spx() are used in MD
          spy(i)= spy00(l)
          spz(i)= spz00(l)
  240     continue
c
        else
          do 243 l= 1,np10
          i= i +1
c
c To avoid always being s=0... 11/23/2010
c  Note: specifying spx(i)= sp2(i)*ranff(0) is not good
          th= pi*ranff(0.)
          ph= pi2*ranff(0.)
          spx(i)= sp2(i)*sin(th)*cos(ph)
          spy(i)= sp2(i)*sin(th)*sin(ph)
          spz(i)= sp2(i)*cos(th)
  243     continue
        end if
  230   continue
c
        if(rank.eq.0) write(06,*) 'np1, np10*cells=',np1,i
        go to 2100
c
      else
        do 250 i= 1,np1
        spx(i)= spx00(i)  ! spx() are used in MD
        spy(i)= spy00(i)
        spz(i)= spz00(i)
c
        aspx(i)= spx00(i)
        aspy(i)= spy00(i)
        aspz(i)= spz00(i)
  250   continue
      end if
c
      kstart= 1
      it= 0
      t= 0
      go to 1000
c
c-----------------------------------------------------------
c* Spin dynamics with dissipation
c    Start from an equilirated state ...
c
c    ds_i                 2*Jij       g*mue_b
c   ------ = s_i x (Sum_j ----- s_j - ------- Be)
c     dt                   hbar        hbar
c
c     f1= J00 *t_unit/hbar        J00= 10 meV
c     g1= g*mue_b*B00 *t_unit/hbar   B00= 100 gauss
c-----------------------------------------------------------
 3000 continue
c                                         ***
      fw= fw00*(1.d0-exp(-t/tau_b)) * sin(omg_b*t)  ! -M*B < 0
      bex= Bapx*fw
      bey= Bapy*fw
      bez= Bapz*fw
c
      iter= 0
  400 iter= iter +1
c
      do 301 i= 1,np1
      if(iter.eq.1) then
        spx1(i)= spx(i)
        spy1(i)= spy(i)
        spz1(i)= spz(i)
      end if
c
      spx0(i)= spx1(i)
      spy0(i)= spy1(i)
      spz0(i)= spz1(i)
  301 continue
c
      do 300 i= i1(rank),i2(rank)
      qsx= 0
      qsy= 0
      qsz= 0
c
      do 320 k= 1,nintS(i)
      j= lintS(k,i)
c
c* Fold back - corrected 11/23/2010
c
      xx= x(i) - x(j)
      yy= y(i) - y(j)
      zz= z(i) - z(j)
      xx= xx -DNINT(xx/xleng)*xleng
      yy= yy -DNINT(yy/yleng)*yleng
      zz= zz -DNINT(zz/zleng)*zleng
c
      rr = sqrt(xx**2 +yy**2 +zz**2)
      J_ki= Jint(k,i)*(r_ij(k,i)/rr)**n_of_J
c
      qsx= qsx +f1*J_ki*(spx1(j) +spx(j))
      qsy= qsy +f1*J_ki*(spy1(j) +spy(j))
      qsz= qsz +f1*J_ki*(spz1(j) +spz(j))
  320 continue !   f1*J_ij * 2s
c
c* Equation is: ds_i/dt= s_i*Q/hbar -(s_i-s_i0)/tau
c                 Q= sum_j (2*Jij*s_j) -g*mue_b*B
c
      qsx= qsx -g1*bex
      qsy= qsy -g1*bey
      qsz= qsz -g1*bez
c
c* RHS
      hh2= dth/tau_diss
      rsx= spx(i) +dth*(spy(i)*qsz -spz(i)*qsy) 
      rsy= spy(i) +dth*(spz(i)*qsx -spx(i)*qsz) 
      rsz= spz(i) +dth*(spx(i)*qsy -spy(i)*qsx) 
     *                      -hh2*(spz(i) -2.d0*aspz(i))
c
c RHS_para
      qq  = qsx**2 +qsy**2 +qsz**2
      rqq = (rsx*qsx +rsy*qsy +rsz*qsz)/qq
c
      dthqq= dth**2 * qq
      rdiv= 1.d0/(1.d0 +dthqq)
c
      spx1(i)= (rsx +dth*(rsy*qsz -rsz*qsy) + dthqq*rqq*qsx)*rdiv
      spy1(i)= (rsy +dth*(rsz*qsx -rsx*qsz) + dthqq*rqq*qsy)*rdiv
      spz1(i)= (rsz +dth*(rsx*qsy -rsy*qsx) + dthqq*rqq*qsz)*rdiv
c
c* Mix value of two iterations
c
      alp= 0.7d0
      spx1(i)= alp*spx1(i) +(1.d0 -alp)*spx0(i)
      spy1(i)= alp*spy1(i) +(1.d0 -alp)*spy0(i)
      spz1(i)= alp*spz1(i) +(1.d0 -alp)*spz0(i)
  300 continue
c
c----------------------
c* Unify the data
c----------------------
      i00= i1(rank) -1
c
      do 360 i= i1(rank),i2(rank)
      wsp1(1,i-i00)= spx1(i) 
      wsp1(2,i-i00)= spy1(i) 
      wsp1(3,i-i00)= spz1(i) 
  360 continue
c
      cnt_send= 3*(i2(rank) -i1(rank) +1)
      call mpi_allgatherv (wsp1,cnt_send,          mpi_real8,
     *                     wsp2,cnt_recv,disp_recv,mpi_real8,
     *                     COMM_LOCAL,ierror)
c
      do 370 i= 1,np1
      spx1(i)= wsp2(1,i)
      spy1(i)= wsp2(2,i)
      spz1(i)= wsp2(3,i)
  370 continue
c
c* Tolerance= 1.d-10 & dt= 2.d-3 (ps)
c   -> iter= 6, error of 0.29% at t= 1000 ps
c
      toler= 1.d-10
      smax= 0
      sav= 0
      notconv= 0
c
      do 390 i= 1,np1
      ss= ((spx1(i) -spx0(i))**2 +(spy1(i) -spy0(i))**2 
     *                           +(spz1(i) -spz0(i))**2)
     *           /(spx0(i)**2 +spy0(i)**2 +spz0(i)**2)
      smax= dmax1(ss,smax)
      sav = sav +ss
c
      if(abs(ss).gt.toler) then
        notconv= notconv +1
      end if
  390 continue
c
      if(notconv.ne.0 .and. iter.lt.itermax) go to 400
c     ------------------------------------------------
c
c* Decrease in spz(i) should be converted to s_perp(i), as sp2(i)= const.
c
      do 500 i= 1,np1
      ss= sp2(i)**2 -spz1(i)**2  ! could be negative by numerics
c
      spmin= 0.05d0*sp2(i)
      if(ss.gt.spmin**2) then
        qq= sqrt(ss/(spx1(i)**2 +spy1(i)**2))
      else
        qq= spmin   ! phase information must be retained
      end if
c
      spx(i)= qq*spx1(i)
      spy(i)= qq*spy1(i)
      spz(i)= spz1(i)
  500 continue
c
c ---------------------------------------------------------------
c     dv_i   tau^2
c   m ---- = ----- Sum_j[ 2nJ_ij(a/r_ij)**(n+1) s_i*s_j] 
c      dt    m*a^2
c                    (e*tau)^2
c                  - --------- K*[x_i -x_i(0)]/a
c                      m*a^3
c
c             f(n)  
c              |
c    v(n-0.5) --- v(n+0.5)
c             x(n) ----- x(n+1)
c
! Move ions
c   mass unit: m_unit is included in fc1, fc2
c ---------------------------------------------------------------
c* Use a sub-time mesh, such that  v*dts < 1
c
      ifdt= 0
      do i= 1,np1+np2
      x0(i)= x(i)
      y0(i)= y(i)
      z0(i)= z(i)
c
      vx0(i)= vx(i)
      vy0(i)= vy(i)
      vz0(i)= vz(i)
      end do
c
c  Slow change ...?
c
      if(mod(it,nt_P3M).eq.1) then
        call P3M_perform (x,y,z,ch,fkx,fky,fkz,np1+np2,iwrt1,
     *                    rank,E_Coulomb_P3M)
      end if
c
 4000 if(ifdt.eq.0) then
        kmax= 1
        dts= dt
      else
        kmax= 2*kmax
        dts= dt/kmax
c
        do i= 1,np1+np2
        x(i)= x0(i)
        y(i)= y0(i)
        z(i)= z0(i)
c
        vx(i)= vx0(i)
        vy(i)= vy0(i)
        vz(i)= vz0(i)
        end do
      end if
c
c
      do 5000 kk= 1,kmax
      call forces (fc1,fc2,fcLJ,alpha,E_sp,E_C_r,E_LJ,
     *             i1,i2,i3,i4,rank,np1,np2,it)
c
c----------------------
c* Unify the force
c----------------------
c 1) Spin force
      i00= i1(rank) -1
c
      do i= i1(rank),i2(rank)
      wsp1(1,i-i00)= fx(i) 
      wsp1(2,i-i00)= fy(i) 
      wsp1(3,i-i00)= fz(i) 
      end do
c
      cnt_send= 3*(i2(rank) -i1(rank) +1)
      call mpi_allgatherv (wsp1,cnt_send,          mpi_real8,
     *                     wsp2,cnt_recv,disp_recv,mpi_real8,
     *                     COMM_LOCAL,ierror)
c
      do i= 1,np1
      fx(i)= wsp2(1,i)
      fy(i)= wsp2(2,i)
      fz(i)= wsp2(3,i)
      end do
c
c 2) Coulomb force
c
      i00= i3(rank) -1
c
      do i= i3(rank),i4(rank)
      wsp1(1,i-i00)= fxC(i) 
      wsp1(2,i-i00)= fyC(i) 
      wsp1(3,i-i00)= fzC(i) 
      end do
c
      cnt_send= 3*(i4(rank) -i3(rank) +1)
      call mpi_allgatherv (wsp1,cnt_send,            mpi_real8,
     *                     wsp2,cnt_recvC,disp_recvC,mpi_real8,
     *                     COMM_LOCAL,ierror)
c
      do i= 1,np1+np2
      fxC(i)= wsp2(1,i)
      fyC(i)= wsp2(2,i)
      fzC(i)= wsp2(3,i)
      end do
c
c 3) Energies
c
      UU1(1)= E_sp
      UU1(2)= E_C_r
      UU1(3)= E_LJ
      call mpi_allreduce (UU1,UU2,3,mpi_real8,mpi_sum,COMM_LOCAL,
     *                    MPIerror)
c
      E_sp = UU2(1)
      E_C_r= UU2(2)
      E_LJ = UU2(3)
c
c
c* Full force should be used from t=0 to avoid atom thermal diffusion
c
      do i= 1,np1+np2
      vx(i)= vx(i) +dts*(fx(i) +fxC(i) +fkx(i))/mass(i)
      vy(i)= vy(i) +dts*(fy(i) +fyC(i) +fky(i))/mass(i)
      vz(i)= vz(i) +dts*(fz(i) +fzC(i) +fkz(i))/mass(i)
c
      x(i)= x(i) +dts*vx(i)
      y(i)= y(i) +dts*vy(i)
      z(i)= z(i) +dts*vz(i)
      end do
c
c  ... check all atoms separately
c
      fdt8= 0
      vdt8= 0
      ifv= 0
c
      do i= 1,np1+np2
      fdt8= dts* dmax1(abs(fx(i)+fxC(i)+fkx(i)),
     *                 abs(fy(i)+fyC(i)+fky(i)),
     *                 abs(fz(i)+fzC(i)+fkz(i)) )/mass(i)
      vdt8= dts* dmax1(abs(vx(i)),abs(vy(i)),abs(vz(i)))
c
      vth= vth0/sqrt(mass(i))
      if(fdt8.gt.0.2*vth .or. vdt8.gt.0.3d0) then
c
        if(ifdt.ge.5 .and. rank.eq.0) then
          if(ifbatch.eq.1) OPEN (unit=wrt,file=praefixc//'.06'//suffix,
     *                 status='old',position='append',form='formatted')
c
           write(wrt,570) it,t,i,ifdt,dts,fdt8/vth,vdt8
           write(wrt,571) x(i),y(i),z(i),vx(i),vy(i),vz(i),
     *                    dts*(fx(i)+fxC(i)+fkx(i))/(mass(i)*vth),
     *                    dts*(fy(i)+fyC(i)+fky(i))/(mass(i)*vth),
     *                    dts*(fz(i)+fzC(i)+fkz(i))/(mass(i)*vth)
  570      format('it,t=',i8,f8.1,'  i,ifdt=',i4,i2,/,
     *            '   dts,fdt8/vth,vdt8=',1pe10.2,0p2f10.5)
  571      format(5x,'x:',1p3e11.3,3x,'vx:',3e11.3,/,5x,'fdt:',3e11.3)
c
           if(ifbatch.eq.1) close (wrt)
        end if
c
        ifv= ifv +1
        if(ifv.ge.1) then
          ifdt= ifdt +1
          go to 4000
        end if
      end if
      end do
 5000 continue 
c
c
c -----------------------------------------------
c* Calculate the system energy (in all steps)
c -----------------------------------------------
c
 5300 continue 
      Ub1= 0
      Um1= 0
c
      do 590 i= i1(rank),i2(rank)
      Ub1= Ub1 + b1 * (bex*spx(i) +bey*spy(i) +bez*spz(i))
c
      do 590 k= 1,nintS(i)
      j= lintS(k,i)
c
c* Fold back
      xx= x(i) - x(j)
      yy= y(i) - y(j)
      zz= z(i) - z(j)
      xx= xx -DNINT(xx/xleng)*xleng
      yy= yy -DNINT(yy/yleng)*yleng
      zz= zz -DNINT(zz/zleng)*zleng
c
      rr = sqrt(xx**2 +yy**2 +zz**2)
      J_ki= Jint(k,i)*(r_ij(k,i)/rr)**n_of_J
c
      Um1= Um1 - 0.5d0*J00* J_ki*(spx(i)*spx(j)
     *                             +spy(i)*spy(j) +spz(i)*spz(j))
  590 continue
c
      UU1(1)= Ub1
      UU1(2)= Um1
      UU1(3)= 0
      call mpi_allreduce (UU1,UU2,3,mpi_real8,mpi_sum,COMM_LOCAL,
     *                    MPIerror)
      Ub= UU2(1)
      Um= UU2(2)
      U1= Ub + Um
c
c ---- diagnosis ----------------- on the major node --------------
c
      if(rank.eq.0) then
c
        if(ifbatch.eq.1) OPEN (unit=wrt,file=praefixc//'.06'//suffix,
     *               status='old',position='append',form='formatted')
c
c Accumulate in every step
        call magnetiz (spx,spy,spz,g,wx1,wy1,wz1,wn1,U1,Uav,wt1,np1,1)
c
c* Energy absorption due to phase lag: 
c       - <M(t)*B(t)/B_0> = - integral M*B/B_0 dt /t 
c
      if(iwrt1.eq.0) then
        ss= 0
        do 735 i= 1,np1
        ss= ss +spx(i)*bex +spy(i)*bey +spz(i)*bez  ! bez= Bz/B00
  735   continue
c
        del_en= del_en - g*ss/np1
c
c ------------------------------
c*  For FT13
c ------------------------------
        OPEN (unit=13,file=praefixc//'.13'//suffix,
     *             status='old',position='append',form='unformatted')
c
        do i= 1,np1
        spx4(i)= spx(i)
        spy4(i)= spy(i)
        spz4(i)= spz(i)
        end do
c
        do i= 1,np1+np2
        ch4(i)= ch(i)
        x4(i)= x(i)
        y4(i)= y(i)
        z4(i)= z(i)
        vx4(i)= vx(i)
        vy4(i)= vy(i)
        vz4(i)= vz(i)
        end do
c
        write(13) t
        write(13) (spx4(i),spy4(i),spz4(i),i=1,np1)
        write(13) (ch4(i),x4(i),y4(i),z4(i),i=1,np1+np2)
        write(13) (vx4(i),vy4(i),vz4(i),i=1,np1+np2)
        close(13)
c
c ------------------------------
c* Make xyz file for DS Viewer
c ------------------------------
c
      OPEN (unit=23,file=praefixc//suffix//'mg.xyz',
     *            status='unknown',position='append',form='formatted')
c                         +++++++
c
      write(23,'(i6)') np1+np2
      write(23,'(a30)') 'All atoms in the entire system'
c
      do i= 1,np1
      write(23,123) 'Fe',x(i),y(i),z(i)
  123 format(a2,3f12.6)
      end do
c
      do i= np1+1,np1+np2
      write(23,123) 'O ',x(i),y(i),z(i)
      end do
c
      close (23)
c
c
c----------------------
c*  History plots
c----------------------
        is= is +1
        if(is.gt.nhs) call rehist(rank)
c
        timeh(is)= t
        do 737 i= 1,np1
        if(spec(i).eq.2) then  ! Fe(2+)
          spinx(is)= spx(i)
          spinz(is)= spz(i)
          spin7(is)= sqrt(spx(i)**2 +spy(i)**2 +spz(i)**2)
          go to 738
        end if
  737   continue
  738   continue
c
        bextx(is)= B00*bex
        bextz(is)= B00*bez
c                                        ! ic= 7: get average
        call magnetiz (spx,spy,spz,g,wx1,wy1,wz1,wn1,U1,Uav,wt1,np1,7)
c
        magx(is)= wx1  ! M= -g*mue_B*s
        magy(is)= wy1  ! U= -M dot B = g*mue_B s*B 
        magz(is)= wz1
c
c
        U_fe(is)= 0
        U_o(is)= 0
        ds1= 0
        ds2= 0
c
        do i= 1,np1
        U_fe(is)= U_fe(is) +0.5*fc3*m_fe*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        ds1= ds1 +sqrt((x(i)-x_0(i))**2 +(y(i)-y_0(i))**2 
     *                                               +(z(i)-z_0(i))**2)
        end do
c
        do i= np1+1,np1+np2
        U_o(is)= U_o(is) +0.5*fc3*m_o*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        ds2= ds2 +sqrt((x(i)-x_0(i))**2 +(y(i)-y_0(i))**2
     *                                               +(z(i)-z_0(i))**2)
        end do
c
        ds_fe(is) = ds1/np1
        ds_o(is)  = ds2/np2
c
        U1 = U1 + U_fe(is) + U_o(is)
c
        Usys8(is)= U1/np1  ! Uav
        Usys(is)=  U1/np1
        Uss(is) =  Um/np1
        Usb(is) =  Ub/np1
c
        U_fe(is)= U_fe(is)/np1
        U_o(is) = U_o(is) /np2
c
        E_sp  = E_sp/np1
        E_C_r = E_C_r/np1
        E_LJ  = E_LJ/np1
c
        t4 = t
        bex4= B00*bex
        bez4= B00*bez
c
c
c--------------------
c* Spin temperature
c--------------------
        do 741 k= 1,3
        sn1(k)= 0
        sx1(k)= 0
        sy1(k)= 0
        sz1(k)= 0
        sx2(k)= 0
        sy2(k)= 0
        sz2(k)= 0
  741   continue
c
        do 743 i= 1,np1
        csx= spx(i)/sqrt(spx(i)**2+spy(i)**2+spz(i)**2)
        csy= spy(i)/sqrt(spx(i)**2+spy(i)**2+spz(i)**2)
        csz= spz(i)/sqrt(spx(i)**2+spy(i)**2+spz(i)**2)
c
        if(spec(i).eq.1) then    ! Fe(3+)
          if(site(i).eq.1) k= 1  !  at (A)
          if(site(i).eq.2) k= 2  !  at (B)
        else
          k= 3   ! Fe(+2) at (B)
        end if
c
        sn1(k)= sn1(k) +1
        sx1(k)= sx1(k) +csx
        sy1(k)= sy1(k) +csy
        sz1(k)= sz1(k) +csz
        sx2(k)= sx2(k) +csx**2
        sy2(k)= sy2(k) +csy**2
        sz2(k)= sz2(k) +csz**2
  743   continue
c
        do 747 k= 1,3
        Tsx(is,k)= sqrt(sx2(k)/sn1(k) -(sx1(k)/sn1(k))**2)
        Tsy(is,k)= sqrt(sy2(k)/sn1(k) -(sy1(k)/sn1(k))**2)
        Tsz(is,k)= sqrt(sz2(k)/sn1(k) -(sz1(k)/sn1(k))**2)
  747   continue
c
c ------------------------------------------------------------
c* av_tsz() is used in the MC/MD step
c
        ss= 0.
        do 748 k= 1,3
        ss= ss +Tsz(is,k)
  748   continue
        av_tsz(is)= ss/3.d0
c
        if(is.eq.5) Tsz0= (av_tsz(3)+av_tsz(4)+av_tsz(5))/3.
c ------------------------------------------------------------
c
        conv(is)= sav/np1
        aitr(is)= iter
        psdt(is)= dtrhs
        sum_mb(is)= del_en
c
        if(is.eq.1) then
          write(wrt,770)
  770     format(//,'# MD run is performed.....',/,
     *      '      it   is       t  itr  Usys      U_fe      U_o     ',
     *      '  conv      F*dts/m_i    v*dt      E_sp     E_C_r    ',
     *      ' E_LJ     magz   dev_x(Fe)   dev_x(O)   cpu(hour)')
        end if
c
        bbw= B00 * bez
        fdt4(is)= fdt8/vth_o
        vdt4(is)= vdt8
        idt4(is)= ifdt
c 
        write(wrt,771) it,is,t,iter,Usys(is),U_fe(is),U_o(is),conv(is),
     *           fdt8,vdt8,E_sp,E_C_r,E_LJ,magz(is),ds_fe(is),ds_o(is),
     *           wtime/3600.d0
  771   format(i8,i5,f9.1,i4,1p6e10.2,1p3e10.2,0pf10.5,2f10.3,f8.2)
c
c                                        ! ic= 0: reset wx-wn
        call magnetiz (spx,spy,spz,g,wx1,wy1,wz1,wn1,U1,Uav,wt1,np1,0)
      end if
c
c
      if(iwrt2.eq.0) then
        OPEN (unit=77,file=praefixc//'.77'//suffix,
     *                status='old',position='append',form='formatted')
        time= t
c
        call plot_spin (x,y,z,spx,spy,spz,spec,site,np1)
        call distr_spin (spx,spy,spz,spec,site,np1)
        close (77)
      end if
c
c
      if(iwrt1.eq.0) then
        OPEN (unit=12,file=praefixc//'.12'//suffix,
     *                            status='unknown',form='unformatted')
c
        write(12) t,xleng,yleng,zleng,rcut,rcutC,Temp,Tcurie,
     *            tmax,dt,cptot
        write(12) x,y,z,vx,vy,vz,ch,mass,ag,ep
        write(12) x_0,y_0,z_0,rintC
        write(12) spx,spy,spz,sp2,spx00,spy00,spz00,r_ij
        write(12) aspx,aspy,aspz
        write(12) Jint,U,spin2,spin3
        write(12) Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,toler,itermax
        write(12) dtwr,dtwr2,fw00,aTsz,Tsz0,av_tsz
        write(12) it,is,iwa,iwb,ir,ir0
        write(12) np1,np2,spec,site,nintS,lintS,nintC,lintC,if_LJ
        write(12) i1,i2,i3,i4,disp_recv,cnt_recv,disp_recvC,cnt_recvC
        write(12) spinx,spinz,spin7,bextx,bextz,magx,magy,magz,
     *            Usys,conv,aitr,psdt,Tfix,Uss,Usb,Tsx,Tsy,Tsz,sum_mb,
     *            U_fe,U_o,fdt4,vdt4,idt4,timeh
        close(12)
      end if
c
      if(ifbatch.eq.1) close (wrt)
      end if
c --------------------------------------- on major nodes --------------
c
        if(istop.eq.1) go to 7000

      GO TO 1000
c
 7000 continue
c
      if(rank.eq.0) then
        OPEN (unit=12,file=praefixc//'.12'//suffix,
     *                           status='unknown',form='unformatted')
c
        write(12) t,xleng,yleng,zleng,rcut,rcutC,Temp,Tcurie,
     *            tmax,dt,cptot
        write(12) x,y,z,vx,vy,vz,ch,mass,ag,ep
        write(12) x_0,y_0,z_0,rintC
        write(12) spx,spy,spz,sp2,spx00,spy00,spz00,r_ij
        write(12) aspx,aspy,aspz
        write(12) Jint,U,spin2,spin3
        write(12) Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,toler,itermax
        write(12) dtwr,dtwr2,fw00,aTsz,Tsz0,av_tsz
        write(12) it,is,iwa,iwb,ir,ir0
        write(12) np1,np2,spec,site,nintS,lintS,nintC,lintC,if_LJ
        write(12) i1,i2,i3,i4,disp_recv,cnt_recv,disp_recvC,cnt_recvC
        write(12) spinx,spinz,spin7,bextx,bextz,magx,magy,magz,
     *            Usys,conv,aitr,psdt,Tfix,Uss,Usb,Tsx,Tsy,Tsz,sum_mb,
     *            U_fe,U_o,fdt4,vdt4,idt4,timeh
        close(12)
c
c
        call date_and_time (cdate,ctime)
c
        if(ifbatch.eq.1) OPEN (unit=wrt,file=praefixc//'.06'//suffix,
     *                status='old',position='append',form='formatted')
c
        write(wrt,*) ' Job is finished  t, cptot=',t,cptot
        write(wrt,*) ' today = ',cdate
        write(wrt,*) ' time  = ',ctime

        write(wrt,*) ' Final: t, tmax=',t,tmax
c
        OPEN (unit=77,file=praefixc//'.77'//suffix,
     *                status='old',position='append',form='formatted')
        call lplots (suffix)
c
        call gclose
c
        close (77)
        if(ifbatch.eq.1) close (wrt)
      end if
c
      return
      end
c
c
c-----------------------------------------------------------------------
      subroutine forces (fc1,fc2,fcLJ,alpha,E_sp,E_C_r,E_LJ,
     *                   i1,i2,i3,i4,rank,np1,np2,it)
c-----------------------------------------------------------------------
      implicit none
c
      include    'param-spinRL.h'
c
      real*8     x,y,z,spx,spy,spz,ch,fx,fy,fz,fxC,fyC,fzC,
     *           Jint,r_ij,rintC,fc1,fc2,fcLJ,alpha,xleng,yleng,zleng,
     *           E_sp,E_C_r,E_LJ,dJdr,r2,r,xx,yy,zz,
     *           ff,pi,sqrtpi,erfc,forceV,rsi,snt,snt2,rlj,ccel
      common/ewald3/ xleng,yleng,zleng
      integer*4  n_of_J,i1(0:num_proc),i2(0:num_proc),i3(0:num_proc),
     *           i4(0:num_proc),rank,nintS,lintS,nintC,lintC,
     *           if_LJ,spec,site,np1,np2,i,j,k,l,it
      common/partcl1/ 
     *           x(np0),y(np0),z(np0),spx(np0),spy(np0),spz(np0),
     *           ch(np0)
      common/partcl2/ 
     *           fx(np0),fy(np0),fz(np0),fxC(np0),fyC(np0),
     *           fzC(np0),Jint(nbx,np0),r_ij(nbx,np0),rintC(nbx,np0)
      common/partcl3/ 
     *           nintS(np0),lintS(nbx,np0),nintC(np0),lintC(nbx,np0),
     *           if_LJ(np0),spec(np0),site(np0)
c
      REAL*8    A1,A2,A3,A4,A5,PP,ar,tt
      PARAMETER ( A1 = 0.254829592d0, A2 = -0.284496736d0 )
      PARAMETER ( A3 = 1.421413741d0, A4 = -1.453152027d0 )
      PARAMETER ( A5 = 1.061405429d0, PP =  0.3275911d0   )
c
      real*8     A,B,C,D,BB
      parameter ( A =  1.8277098d3,     ! e^2/A
     *            B =  4.9254591d0,     ! 1/A
     *            C = -2.1362173d0,     ! e^2/A^5
     *            D =  7.4679860d1  )   ! e^2/A^11
c*---------------------------------------------------------------
c
      pi= 4*atan(1.d0)
      sqrtpi= sqrt(pi)
c
      do i= 1,np1+np2
      fx(i)= 0   ! note: not defined for i > np1
      fy(i)= 0
      fz(i)= 0
      end do
c
      do i= i3(rank),i4(rank)
      fxC(i)= 0
      fyC(i)= 0
      fzC(i)= 0
      end do
c
      E_sp = 0
      E_LJ = 0
      E_C_r= 0
c
c Spin force
c
      do i= i1(rank),i2(rank)
      do k= 1,nintS(i)
      j= lintS(k,i)
c
! 2*0.5 is omitted in dJdr
c     fc1= J00*(t_unit/a_unit)**2/m_unit          ! spin force
c     fc2= (e_unit*t_unit)**2/(m_unit*a_unit**3)  ! Coulomb force
c
c* Fold back
      xx= x(i) - x(j)
      yy= y(i) - y(j)
      zz= z(i) - z(j)
      xx= xx -DNINT(xx/xleng)*xleng
      yy= yy -DNINT(yy/yleng)*yleng
      zz= zz -DNINT(zz/zleng)*zleng
c
      r = sqrt(xx**2 +yy**2 +zz**2)
      dJdr= -fc1 *Jint(k,i)*(r_ij(k,i)/r**2)*
     *            (spx(i)*spx(j) +spy(i)*spy(j) +spz(i)*spz(j))
c
      fx(i)= fx(i) - dJdr*xx/r
      fy(i)= fy(i) - dJdr*yy/r
      fz(i)= fz(i) - dJdr*zz/r
c
      E_sp= E_sp +dJdr
      end do
      end do
c
c
c* Coulomb and LJ forces
c
      do i= i3(rank),i4(rank)
      do k= 1,nintC(i)
      j= lintC(k,i)
c
      xx= x(i) - x(j)
      yy= y(i) - y(j)
      zz= z(i) - z(j)
      xx= xx -DNINT(xx/xleng)*xleng
      yy= yy -DNINT(yy/yleng)*yleng
      zz= zz -DNINT(zz/zleng)*zleng
c
      r2= xx**2 +yy**2 +zz**2
      r = sqrt(r2)
c
      rsi = 1/r2
      snt = rsi*rsi*rsi
c
c*  Lennard-Jones force only for the Fe-O pair
c     ref. J.Rustad, B.Hay, and J.Halley, J.Chem.Phys. 102, 427 (1995).
c             ***
c
      if(spec(i).eq.0) then         ! O-
        if(site(j).eq.1) go to 300  !  -Fe(A)
        if(site(j).eq.2) go to 400  !  -Fe(B) +3, +2
      else                            ! Fe-
        if(spec(j).eq.0) then         !   -O
          if(site(i).eq.1) go to 300  ! Fe on A-site
          if(site(i).eq.2) go to 400  ! Fe on B-site
        end if
      end if
      go to 500
c
  300 continue
      E_LJ = E_LJ +fc2*( A*exp(-B*r) +(C +D*snt)*snt )
      ccel =  fc2*( A*B*exp(-B*r)/r +(6*C +12*D*snt)*snt/r2 )
      go to 600
c
  400 continue
      if(.true.) go to 300
c     +++++++++ 
c
      if(spec(i).eq.1 .or. spec(j).eq.1) then
        BB = 0.85d0*B  ! +3
      else
        BB = 0.9d0*B   ! +2
      end if
c
      E_LJ = E_LJ +fc2*( A*exp(-BB*r) +(C +D*snt)*snt )
      ccel =  fc2*( A*BB*exp(-BB*r)/r +(6*C +12*D*snt)*snt/r2 )
      go to 600
c
c  Fe-Fe, or O-O pair
  500 continue
c
      ccel= 0
      if(if_LJ(i).eq.1) then   ! limit the pair to r(t=0) < rcutLJ
c        ++++++++
c
        snt2= 2**(1.d0/6.d0)
        rlj = r/(rintC(k,i)/snt2)   ! minimum at rintC(k,i) for this pair
c                *****
        rsi = 1/dmax1(rlj**2, 0.81d0)
        snt = rsi*rsi*rsi
c
        E_LJ = E_LJ +fc2*(D/12.d0)*snt*(snt -1.d0)
        ccel =       fc2*D*snt*(snt -0.5d0)/r2
      end if
c
  600 continue
c
c*  Coulomb force in Ewald sum
c
      ar   = alpha*r
      tt   = 1.d0 / ( 1.d0 + PP * ar )
      erfc = tt*(A1+tt*(A2+tt*(A3+tt*(A4+tt*A5))))
c
      E_C_r  = E_C_r + fc2*ch(i)*ch(j)* erfc*dexp(-ar**2)/r
      forceV = fc2*ch(i)*ch(j)*(erfc/r +2*alpha/sqrtpi)*dexp(-ar**2)/r2
c        1         2         3         4         5         6         71
c
      fxC(i)= fxC(i) +(ccel +forceV)*xx 
      fyC(i)= fyC(i) +(ccel +forceV)*yy
      fzC(i)= fzC(i) +(ccel +forceV)*zz
      end do
      end do
c
      return
      end
c
c
c-----------------------------------------------------------------
      subroutine init (x,y,z,ch,spx,spy,spz,sp2,spec,site,
     *                 xleng0,yleng0,zleng0,np1,np2,np10,
     *                 i1x,i2x,i1y,i2y,i1z,i2z,rank,size,suffix)
c-----------------------------------------------------------------
cO 2.057265000 2.057265000 2.057265000
      implicit none
c
      include    'mpif.h'
      include    'param-spinRL.h'
c
      real*8     xo(100),yo(100),zo(100),cho(100)
      real*8     x(np0),y(np0),z(np0),ch(np0),spx(np0),spy(np0),
     *           spz(np0),sp2(np0),spin2,spin3,Jaa,Jbb,Jab,
     *           Bapx,Bapy,Bapz,tau_b,tau_diss,Temp,Tcurie,
     *           xleng,yleng,zleng,rcut,pi,pi2,th,ph,xx,yy,zz,
     *           t,dt,tmax,cptot,xleng0,yleng0,zleng0,
     *           a(3),b(3),c(3),qfe2,qfe3,qo
      real*8     eps,ranff
      integer*4  spec(np0),site(np0),np1,np2,rank,size,np10,np20,
     *           mx,my,mz,i,l,j,lo,ia,ja,ka,ic,jc,kc,i1,i2,j1,j2,k1,k2,
     *           nfe2,nfe3,ifbatch,wrt,ierror,
     *           i1x,i2x,i1y,i2y,i1z,i2z
      character  suffix*1,char*2,text1*125
c
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,
     *              tau_diss,Temp,Tcurie
      common/parm3/ t,rcut,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
      common/parm4/ mx,my,mz
      common/print6/ ifbatch,wrt
c
c* Three axes of the crystal
      a(1)= xleng0
      a(2)= 0.000000000 
      a(3)= 0.000000000 
      b(1)= 0.000000000 
      b(2)= yleng0
      b(3)= 0.000000000 
      c(1)= 0.000000000 
      c(2)= 0.000000000 
      c(3)= zleng0
c
      if(rank.eq.0) then
        OPEN (unit=17,file='magnetite8.xyz',
     *                                status='old',form='formatted')
c
        read(17,'(i6)') np10 
        read(17,'(a125)') text1
      end if
      call mpi_bcast ( np10,  1,mpi_integer,  0,MPI_COMM_WORLD,ierror)
c     call mpi_bcast (text1,125,mpi_character,0,MPI_COMM_WORLD,ierror)
c
      pi= 4*datan(1.d0)
      pi2= 2*pi
c
      spin2= 2.0d0  ! Fe 2+
      spin3= 2.5d0  ! Fe 3+
c
      qfe2= 2.0d0
      qfe3= 3.0d0
      qo = -2.0d0
c
      i= 0
      j= 0
      l= 0
c
c  magnet3.xyz
c  FT17: O(1)...O(32) B(1)...B(16) A(1)...A(8) 
c
  100 l= l +1
      if(l.gt.np10) go to 700
c
c* Balance charges on cubic sub-lattice
c
      if(l.gt.32 .and. mod(l,4).eq.1) then
        nfe2= 0
        nfe3= 0
      end if
c
      if(rank.eq.0) then
        read(17,101) char,xx,yy,zz
  101   format(a2,f12.0,f12.0,f12.0)
      end if
c
      call mpi_bcast (char,2,mpi_character,     0,MPI_COMM_WORLD,ierror)
      call mpi_bcast (xx,1,mpi_double_precision,0,MPI_COMM_WORLD,ierror)
      call mpi_bcast (yy,1,mpi_double_precision,0,MPI_COMM_WORLD,ierror)
      call mpi_bcast (zz,1,mpi_double_precision,0,MPI_COMM_WORLD,ierror)
c
! Oxygen
      if(l.le.32) then  
        j= j +1
        xo(j)= xx
        yo(j)= yy
        zo(j)= zz
        cho(j)= qo
        spec(j)= 0     ! oxygen
        go to 100
      end if
c
! Fe
      i= i+1
      x(i)= xx
      y(i)= yy
      z(i)= zz
c
      if(l.ge.49 .and. l.le.56) then
        ch(i)= qfe3
        sp2(i)= spin3
        spec(i)= 1
        site(i)= 1      ! A site (8)
      else
        eps= ranff(0.)
        if(eps.gt.0.5d0) then
          if(nfe3.lt.2) go to 110  ! ...2 Fe(3+) are allowed on sub lattice
          go to 111
        else
          if(nfe2.lt.2) go to 111  ! ...2 Fe(2+) are allowed
          go to 110
        end if
c
  110   nfe3= nfe3 +1
        ch(i)= qfe3
        sp2(i)= spin3
        spec(i)= 1       ! 3+
        site(i)= 2       ! B site (8+8)
        go to 130
c
  111   nfe2= nfe2 +1
        ch(i)= qfe2
        sp2(i)= spin2
        spec(i)= 2       ! 2+
        site(i)= 2       ! B site (8+8)
c
  130   continue
      end if
c  
      th= pi*ranff(0.)
      ph= pi2*ranff(0.)
      spx(i)= sp2(i)*sin(th)*cos(ph)
      spy(i)= sp2(i)*sin(th)*sin(ph)
      spz(i)= sp2(i)*cos(th)
      go to 100
c
  700 np10= i
      np20= j

      if(rank.eq.0) then
        close (17)
c
        write(wrt,703) np10,np20,nfe2,nfe3+8
  703   format(/,' Fundamental domain: #Fe:#O, Fe(2+):Fe(3+)=',2i4,2i4)
      end if
c
c------------------------
c* Expand the system
c------------------------
c  (mx,my,mz) the number of domains - must be odd numbers
c
      l= np10  ! ... Fe
c
      xleng= mx*xleng0
      yleng= my*yleng0
      zleng= mz*zleng0
c
      if(mx.eq.1) then
        i1= 0
        i2= 0
      else
        if(mod(mx,2).eq.1) then
          i1= -mx/2
          i2=  mx/2
        else
          i1= -mx/2 +1
          i2=  mx/2
        end if
      end if
c
      if(my.eq.1) then
        j1= 0
        j2= 0
      else
        if(mod(my,2).eq.1) then
          j1= -my/2
          j2=  my/2
        else
          j1= -my/2 +1
          j2=  my/2
        end if
      end if
c
      if(mz.eq.1) then
        k1= 0
        k2= 0
      else
        if(mod(mz,2).eq.1) then
          k1= -mz/2
          k2=  mz/2
        else
          k1= -mz/2 +1
          k2=  mz/2
        end if
      end if
c
      i1x= i1
      i2x= i2
      i1y= j1
      i2y= j2
      i1z= k1
      i2z= k2
c
      if(mx.eq.1) go to 301
c     +++++++++++++++++++++  i1=i2=0 N/G on Infiniband 
c
      do 300 ia= i1,i2
      do 300 ja= j1,j2
      do 300 ka= k1,k2
      if(ia.eq.0 .and. ja.eq.0 .and. ka.eq.0) go to 300
c        ++++++++ fundamental domain - do not change !!
c
! Fe
      do 400 i= 1,np10
      l = l +1
      x(l)= x(i) +ia*a(1) +ja*b(1) +ka*c(1)
      y(l)= y(i) +ia*a(2) +ja*b(2) +ka*c(2)
      z(l)= z(i) +ia*a(3) +ja*b(3) +ka*c(3)
c
c  Fe: B(1)...B(16) A(17)...A(24) 
      if(i.le.16 .and. mod(i,4).eq.1) then
        nfe2= 0
        nfe3= 0
      end if
c
! A site
      if(i.gt.16) then  
        ch(l) = ch(i)
        sp2(l) = sp2(i)
        spec(l)= spec(i)
        site(l)= site(i)
        go to 430
      else
c
! B site
c* Shuffle Fe(2+)/Fe(3+)= 8/8
c
        if(ranff(0.).gt.0.5d0) then
          if(nfe3.lt.2) go to 410  ! accomodate 2 on sub-lattice
          go to 411
        else
          if(nfe2.lt.2) go to 411
          go to 410
        end if
      end if
c
  410 nfe3= nfe3 +1
      ch(l) = qfe3
      sp2(l)= spin3
      spec(l)= 1
      site(l)= 2       ! B site (16)
      go to 430
c
  411 nfe2= nfe2 +1
      ch(l) = qfe2
      sp2(l)= spin2
      spec(l)= 2
      site(l)= 2       ! B site (16)
c
  430 continue
      th= pi*ranff(0.)
      ph= pi2*ranff(0.)
      spx(l)= sp2(i)*sin(th)*cos(ph)
      spy(l)= sp2(i)*sin(th)*sin(ph)
      spz(l)= sp2(i)*cos(th)
  400 continue
  300 continue
c
  301 continue
c --------------
      np1= l
c --------------
c
! Oxygen
c
      j= 0
c
      do 500 i= 1,np20
      l= l +1
      j= j +1
      x(l)= xo(i)
      y(l)= yo(i)
      z(l)= zo(i)
      ch(l)= qo
      spec(l)= 0
  500 continue
c
      if(mx.eq.1) go to 601
c     +++++++++++++++++++++  i1=i2=0 N/G on Infiniband 
      do 600 ia= i1,i2
      do 600 ja= j1,j2
      do 600 ka= k1,k2
      if(ia.eq.0 .and. ja.eq.0 .and. ka.eq.0) go to 600
c
      do 630 i= 1,np20
      l= l +1
      j= j +1
      x(l)= xo(i) +ia*a(1) +ja*b(1) +ka*c(1)
      y(l)= yo(i) +ia*a(2) +ja*b(2) +ka*c(2)
      z(l)= zo(i) +ia*a(3) +ja*b(3) +ka*c(3)
      ch(l)= qo
      spec(l)= 0
  630 continue
c
  600 continue
  601 continue
c
c --------------
      np2= j
c --------------
c
c  File i/o is allowed only by the master node
c
      if(rank.eq.0) then
c
        write(wrt,707) size,np1,np2
  707   format(' init: num of procs=',i3,'  np(Fe),np(O)=',2i6,/)
c
        close (17)
c +++++
        OPEN (unit=18,file='mag3.xyz',status='unknown',
     *                                      form='formatted')
c
        write(18,'(i6)') np1,np2
        write(18,'(a80)') 'Magnetite in 3d'
c
        do 105 i= 1,np1
        j= np1 +i
        write(18,107) 'Fe ',i,x(i),y(i),z(i),ch(i),
     *                ' O ',j,x(j),y(j),z(j),ch(j)
  107   format(a3,'i=',i4,3f10.5,' q=',f5.1,2x,
     *         a3,'i=',i4,3f10.5,' q=',f5.1)
  105   continue
c
        do 115 i= np1+1,np2
        j= np1 +i
        write(18,108) ' O ',j,x(j),y(j),z(j),ch(j)
  108   format(49x,a3,'i=',i4,3f10.5,' q=',f5.1)
  115   continue
c
        close (18)
c
c* Make xyz file for DS Viewer
c
        OPEN (unit=23,file=praefixc//suffix//'mg.xyz',
     *                           status='unknown',form='formatted')
c
        write(23,'(i6)') np1+np2
        write(23,'(a30)') 'All atoms in the entire system'
c
        do i= 1,np1
        write(23,123) 'Fe',x(i),y(i),z(i)
  123   format(a2,3f12.6)
        end do
c
        do i= np1+1,np1+np2
        write(23,123) 'O ',x(i),y(i),z(i)
        end do
c
        close (23)
c +++++
      end if
c
      return
      end
c
c
CCCCCCCCCCCCCCC   P3M (FORTRAN 77)  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c                                 26.10.99 
c                                     Motohiko Tanaka, Christian Holm
c  Caling sequences:
c
c     call  P3M_init (length,alpha0,IP)
c     call  P3M_perform (coox,cooy,cooz,Q,FX,FY,FZ,E_Coulomb_P3M)
c
c/*---------------------------------------------------------------------
c SUBUNIT:  p3m_v2.f  (FORTRAN 90)
C 
C       Version   20.10.1999 (CH) 
C       corrected 26.10.1999 (MT)
c                 23.06.2000 (MT)
c
c VERSION:  20 January 1999
c AUTHOR:   Markus Deserno
c
c    Brillouin is now a parameter
c    MaxInterpol --> MINTPOL
c    floor --> Defined: must be stepwise at x= 0.
c    dround --> = DNINT (Round off to nearest integer). 
c/*---------------------------------------------------------------------
      subroutine P3M_perform (x,y,z,ch,fkx,fky,fkz,np,iwrt1,
     *                        rank,E_Coulomb_P3M)
c/*---------------------------------------------------------------------
      IMPLICIT    NONE
c
      include     'mpif.h'
      include     'fftw3.f'
c
      include     'param-spinRL.h'
      integer*4   np,iwrt1,rank,ierror
c
      real*8      x(0:np0-1),y(0:np0-1),z(0:np0-1),ch(0:np0-1),
     *            fkx(0:np0-1),fky(0:np0-1),fkz(0:np0-1),
     *            ecl,E_Coulomb_P3M,wupi,pi
c
      integer*4   xarg(0:np0-1),yarg(0:np0-1),zarg(0:np0-1), 
     *            i,j,k,m,MESHMASKx,MESHMASKy,MESHMASKz,
     *            Gi0,Gi1,Gi2,xpos,ypos,zpos,
     *            assignshiftx,assignshifty,assignshiftz,
     *            QZahl,m0
      real*8      d1,Hix,Hiy,Hiz,MI2,modadd1,modadd2,T1,T2,T3,
     *            sum_q_2, sum_q2
c     
c-----------
      real*8         alpha,fc2,xleng,yleng,zleng
      real*8         meshiftx,meshifty,meshiftz,Dnx,Dny,Dnz
      integer*4      PP
      common/ewald1/ alpha,fc2
      common/ewald2/ PP
      common/ewald3/ xleng,yleng,zleng
c
      real*8         cooP,QP,QL,Ghat,intCAF
      integer*4      global,G
c-----------
      common/coordi/ cooP(0:np0-1,0:2),QP(0:np0-1)
      common/pindex/ global(0:np0-1),G(0:np0-1,0:2)
      common/prmesh/ QL(0:P_max**3-1,0:np0-1)
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy),
     *               meshiftz(0:meshz),Dnx(0:meshx),Dny(0:meshy),
     *               Dnz(0:meshz)
c
      common/influf/ Ghat(0:meshx,0:meshy,0:meshz)
      common/intcaf/ intCAF(0:P_max,0:2*MINTPOL+1)
c
      complex*16         qq_c,phi_x_c,phi_y_c,phi_z_c
      common/fft001/     qq_c(0:meshx/2, 0:meshy-1,0:meshz-1),  ! FFTW-style
     *                phi_x_c(0:meshx/2, 0:meshy-1,0:meshz-1),  !  arrays
     *                phi_y_c(0:meshx/2, 0:meshy-1,0:meshz-1),
     *                phi_z_c(0:meshx/2, 0:meshy-1,0:meshz-1)
c
      real*8             qq,phi_x,phi_y,phi_z
      common/fft002/     qq(0:meshx-1,0:meshy-1,0:meshz-1),  
     *                phi_x(0:meshx-1,0:meshy-1,0:meshz-1),  
     *                phi_y(0:meshx-1,0:meshy-1,0:meshz-1),  
     *                phi_z(0:meshx-1,0:meshy-1,0:meshz-1)  
c-----------
c
c -------------------------
c*  FFTW - FORTRAN header
c -------------------------
c       EM64T: 64-bit pointer
c                 **** ****
      integer*8   plan,pinv
c
      logical     first
      common/fftsav/ plan,pinv,first
      data        first /.true./
c
      real*8      fft_scale2
      complex*16  ei,eigq
c
c ---------------------------
c*  Prepare for FFTW calls
c ---------------------------
c
      if(first) then
        first= .false.
c
c       call dfftw_init_threads (ierror)
c       call dfftw_plan_with_nthreads (nthreads)
c
        call dfftw_plan_dft_r2c_3d (plan,meshx,meshy,meshz,qq,qq_c,
     *                              FFTW_ESTIMATE)
c
        call dfftw_plan_dft_c2r_3d (pinv,meshx,meshy,meshz,qq_c,qq,
     *                              FFTW_ESTIMATE)
      end if
c
c
      ei = dcmplx(0.d0, 1.d0)
c --------------------------------
      pi = 4.d0*DATAN(1.d0)
      MESHMASKx = meshx-1     
      MESHMASKy = meshy-1     
      MESHMASKz = meshz-1     
c
      Hix = meshx / xleng
      Hiy = meshy / yleng
      Hiz = meshz / zleng
c
      MI2 = 2.d0*dfloat(MINTPOL)
      assignshiftx = meshx -(PP -1)/2
      assignshifty = meshy -(PP -1)/2
      assignshiftz = meshz -(PP -1)/2
c
      QZahl = 0
      sum_q_2 = 0.d0
      sum_q2  = 0.d0
c
      do 100 i= 0,np-1
      if (dabs(ch(i)) .gt. 1.d-5) then 
        cooP(QZahl, 0) = x(i) - DNINT(x(i)/xleng -0.5d0)*xleng 
        cooP(QZahl, 1) = y(i) - DNINT(y(i)/yleng -0.5d0)*yleng
        cooP(QZahl, 2) = z(i) - DNINT(z(i)/zleng -0.5d0)*zleng
c
        QP(QZahl) = ch(i)
        sum_q_2 = sum_q_2 + QP(QZahl)
        sum_q2  = sum_q2  + QP(QZahl)**2
c
        global(QZahl) = i
        QZahl= QZahl + 1
      end if
c
      fkx(i)= 0
      fky(i)= 0
      fkz(i)= 0
  100 continue
c
      sum_q_2 = sum_q_2 **2
c
      do 200 i= 0,meshx-1 
      do 200 j= 0,meshy-1 
      do 200 k= 0,meshz-1 
      qq(i,j,k) = 0.d0
  200 continue
c
      if(PP.eq.3) then
         modadd1 = 0.0d0
         modadd2 = 0.5d0
      else
         write(6,*) "Error in function 'P3M_perform':"
         write(6,*) "Charge assignment order P=",PP," unknown."
         write(6,*) "Program terminated."
         call exit(1)
      end if
c
c
      do 300 i= 0, QZahl-1
c*
      d1  = cooP(i,0)*Hix + modadd1
      Gi0 = int(d1 + modadd2) + assignshiftx
      G(i,0) = Gi0 
      xarg(i) = int( (d1 - dnint(d1) + 0.5d0)*MI2 )
c      
      d1  = cooP(i,1)*Hiy + modadd1 
      Gi1 = int(d1 + modadd2) + assignshifty
      G(i,1) = Gi1
      yarg(i) = int( (d1 - dnint(d1) + 0.5d0)*MI2 )
c      
      d1  = cooP(i,2)*Hiz + modadd1 
      Gi2 = int(d1 + modadd2) + assignshiftz
      G(i,2) = Gi2 
      zarg(i) = int( (d1 - dnint(d1) + 0.5d0)*MI2 )
c
      m0= -1
c*
      do 330 j = 0, PP-1
      do 330 k = 0, PP-1
      do 330 m = 0, PP-1
      T1 = QP(i) * intCAF(j,xarg(i))
      T2 = T1 *    intCAF(k,yarg(i))
      T3 = T2 *    intCAF(m,zarg(i))
c    
      m0= m0 + 1
      QL(m0,i) = T3          ! assignment factor.
c      
      xpos = IAND( (G(i,0) + j), MESHMASKx)
      ypos = IAND( (G(i,1) + k), MESHMASKy)
      zpos = IAND( (G(i,2) + m), MESHMASKz)
c    
      qq(xpos,ypos,zpos) = qq(xpos,ypos,zpos) + QL(m0,i)
  330 continue
  300 continue
c
c
      fft_scale2= dfloat(meshx*meshy*meshz)
      call dfftw_execute_dft_r2c (plan,qq,qq_c)
c
c*
      ecl = 0.d0
c
      do 700  i= 0, meshx/2
      do 700  j= 0, meshy-1
      do 700  k= 0, meshz-1
      ecl = ecl + Ghat(i,j,k)*cdabs(qq_c(i,j,k))**2/fft_scale2 
  700 continue
c
c-----------------------------------------------------------------
c*  Prefactor is for FFT3R.    
c      checked by Madelung constant                 10/06/2000
c-----------------------------------------------------------------
c
      E_Coulomb_P3M = fc2 * ecl* (1.d0/meshx) *xleng/(4.d0*pi)
c
c     wupi = dsqrt(pi)
c     E_Coulomb_P3M = E_Coulomb_P3M - 
c    *         L *prefactor *( sum_q2 *alpha /wupi
c    *                       + sum_q_2 *pi /(2.d0*L**3*alpha**2) )
c
c
      do 800 j= 0,meshy-1
      do 800 k= 0,meshz-1
      phi_x_c(0,j,k) = 0.d0
      phi_y_c(0,j,k) = 0.d0
      phi_z_c(0,j,k) = 0.d0
  800 continue
c
      do 810 i= 0,meshx/2
      do 810 j= 0,meshy-1
      do 810 k= 0,meshz-1
      eigq = ei*Ghat(i,j,k)*qq_c(i,j,k)/fft_scale2
c
      phi_x_c(i,j,k) = eigq*Dnx(i)
      phi_y_c(i,j,k) = eigq*Dny(j)
      phi_z_c(i,j,k) = eigq*Dnz(k)
  810 continue
c    
c*  For diagnosis.
c    
      do 860 i= 0,meshx/2
      do 860 j= 0,meshy-1
      do 860 k= 0,meshz-1
      qq_c(i,j,k) = Ghat(i,j,k)*qq_c(i,j,k)/fft_scale2
  860 continue
c
      if(iwrt1.eq.0) then
        call dfftw_execute_dft_c2r (plan,qq_c,qq)
      end if
c
      call dfftw_execute_dft_c2r (plan,phi_x_c,phi_x)
      call dfftw_execute_dft_c2r (plan,phi_y_c,phi_y)
      call dfftw_execute_dft_c2r (plan,phi_z_c,phi_z)
c
c
      do 900 i = 0, QZahl-1
c
      m0= -1
      do 930 j = 0, PP-1
      do 930 k = 0, PP-1
      do 930 m = 0, PP-1
      xpos = IAND( (G(i,0) + j), MESHMASKx)
      ypos = IAND( (G(i,1) + k), MESHMASKy)
      zpos = IAND( (G(i,2) + m), MESHMASKz)
c    
      m0 = m0 +1
      d1 = fc2 * QL(m0,i)
c    
      fkx(global(i)) = fkx(global(i)) - d1*phi_x(xpos,ypos,zpos)
      fky(global(i)) = fky(global(i)) - d1*phi_y(xpos,ypos,zpos)
      fkz(global(i)) = fkz(global(i)) - d1*phi_z(xpos,ypos,zpos)
  930 continue
  900 continue
c
      return
      end
c
c
c/*---------------------------------------------------------------------
      subroutine perform_aliasing_sums (NX,NY,NZ,nominatorX,nominatorY, 
     *                                  nominatorZ,denominator)
c/*---------------------------------------------------------------------
      IMPLICIT    NONE
c
      include    'param-spinRL.h'
      integer*4   Brillouin
      parameter  (Brillouin=1)
c
      integer*4      PP
      real*8         alpha,fc2,xleng,yleng,zleng
      real*8         meshiftx,meshifty,meshiftz,Dnx,Dny,Dnz,
     *               dmeshx,dmeshy,dmeshz,pi
      common/ewald1/ alpha,fc2
      common/ewald2/ PP
      common/ewald3/ xleng,yleng,zleng
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy),
     *               meshiftz(0:meshz),Dnx(0:meshx),Dny(0:meshy),
     *               Dnz(0:meshz)
c
      integer*4   NX,NY,NZ,MX,MY,MZ
      real*8      nominatorX,nominatorY,nominatorZ,denominator,
     *            S1,S2,S3,fak1x,fak1y,fak1z,fak2,fak3,
     *            NMX,NMY,NMZ,NM2,expo,exponent_limit,sinc
c
      exponent_limit = 30.d0
c
      dmeshx= dfloat(meshx)
      dmeshy= dfloat(meshy)
      dmeshz= dfloat(meshz)
c
      pi = 4.d0*DATAN(1.d0)
c
c     fak2 = (pi/alpha)**2/Leng**2
      fak2 = (pi/alpha)**2
c
      nominatorX = 0.d0
      nominatorY = 0.d0
      nominatorZ = 0.d0
      denominator = 0.d0
c
      do 100 MX = -Brillouin, Brillouin
      NMX = meshiftx(NX) + dmeshx*MX
      S1  = sinc(NMX/dmeshx)**(2*PP)
c
        do 200 MY = -Brillouin, Brillouin
        NMY = meshifty(NY) + dmeshy*MY
        S2  = S1* sinc(NMY/dmeshy)**(2*PP)
c
          do 300 MZ = -Brillouin, Brillouin
          NMZ = meshiftz(NZ) + dmeshz*MZ
          S3  = S2* sinc(NMZ/dmeshz)**(2*PP)
c
          denominator = denominator + S3
          NM2 = NMX**2 + NMY**2 + NMZ**2
c
          expo= fak2*((NMX/xleng)**2 +(NMY/yleng)**2 +(NMZ/zleng)**2)
          if (expo .lt. exponent_limit) then
              fak3 =  S3* dexp(-expo)/NM2 
          else
              fak3 = 0.d0
          end if
c
          nominatorX = nominatorX + fak3 * NMX
          nominatorY = nominatorY + fak3 * NMY
          nominatorZ = nominatorZ + fak3 * NMZ
  300     continue
  200   continue
  100 continue
c
      return
      end
c
c
c/*---------------------------------------------------------------------
      subroutine  calculate_differential_operator (rank,wrt)
c/*---------------------------------------------------------------------
      IMPLICIT    NONE
c
      include     'param-spinRL.h'
      integer*4   rank,wrt,i,j,k
c-----------
      real*8         meshiftx,meshifty,meshiftz,Dnx,Dny,Dnz,
     *               dmeshx,dmeshy,dmeshz
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy),
     *               meshiftz(0:meshz),Dnx(0:meshx),Dny(0:meshy),
     *               Dnz(0:meshz)
c                                       ******  defined here.
      if(rank.eq.0) then
        write(wrt,*) ' - calculating calculate_differential_operator'
      end if
c
      dmeshx= dfloat(meshx)
      dmeshy= dfloat(meshy)
      dmeshz= dfloat(meshz)
c
      do i= 0,meshx-1
      Dnx(i) = dfloat(i) - dnint(dfloat(i)/dmeshx)*dmeshx
      end do
c
      do j= 0,meshy-1
      Dny(j) = dfloat(j) - dnint(dfloat(j)/dmeshy)*dmeshy
      end do
c
      do k= 0,meshz-1
      Dnz(k) = dfloat(k) - dnint(dfloat(k)/dmeshz)*dmeshz
      end do
c
      Dnx(meshx/2) = 0.d0
      Dny(meshy/2) = 0.d0
      Dnz(meshz/2) = 0.d0
c
      return
      end
c
c
c/*---------------------------------------------------------------------
      subroutine  calculate_influence_function (rank,wrt)
c/*---------------------------------------------------------------------
      IMPLICIT    NONE
c
      include     'param-spinRL.h'
      integer*4    rank,wrt
c-----------
      real*8         alpha,fc2,Ghat,xleng,yleng,zleng,Leng2
      real*8         meshiftx,meshifty,meshiftz,Dnx,Dny,Dnz
      integer*4      PP
      common/ewald1/ alpha,fc2
      common/ewald2/ PP
      common/ewald3/ xleng,yleng,zleng
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy),
     *               meshiftz(0:meshz),Dnx(0:meshx),Dny(0:meshy),
     *               Dnz(0:meshz)
      common/influf/ Ghat(0:meshx,0:meshy,0:meshz)
c                    **************************  defined here.
c
      integer*4   NX,NY,NZ
      real*8      DDnx,DDny,DDnz,Dn2,fak1,fak2,fak3,
     *            nominatorX,nominatorY,nominatorZ,denominator
c
      if(rank.eq.0) then
        write(wrt,*) ' - calculating influence function with parameters'
        write(wrt,601) meshx,meshy,meshz,PP,alpha,xleng,yleng,zleng
  601   format('   meshx, meshy, meshz=',3i5,/,
     *         '   P=',i5,/,'   alpha=',d20.12,/,
     *         ',  xleng, yleng, zleng=',3d20.12,/)
      end if
c 
c     fak1= 2.d0*dfloat(mesh*mesh*mesh)/Leng**2
c
c         write(11,*) 'P3M - Aliasing sums for different NX,NY,NZ:'
c  
      do 100 NX = 0, meshx-1
       do 200 NY = 0, meshy-1
        do 300 NZ = 0, meshz-1
        if ( (NX.eq.0).and.(NY.eq.0).and.(NZ.eq.0)) then
           Ghat(NX,NY,NZ)= 0.d0 
        else
           call perform_aliasing_sums (NX,NY,NZ,nominatorX,nominatorY,
     *                                 nominatorZ,denominator)
           DDnx = Dnx(NX)
           DDny = Dny(NY)
           DDnz = Dnz(NZ)
c
           Dn2 = DDnx**2 + DDny**2 + DDnz**2
c  
           if (Dn2 .gt. 1.d-7) then
             Ghat(NX,NY,NZ) = 2.d0*( 
     *                     meshx*DDnx*nominatorX/xleng**2 
     *                   + meshy*DDny*nominatorY/yleng**2 
     *                   + meshz*DDnz*nominatorZ/zleng**2 )
     *                   /( Dn2 * denominator**2 )
           else 
             Ghat(NX,NY,NZ) = 0.d0
           end if
        end if
c
  300   continue
  200  continue
  100 continue
c
      return
      end
c
c
c/*---------------------------------------------------------------------
      subroutine interpol_charge_assign_function (rank,wrt)
c/*---------------------------------------------------------------------
c*  P = 3 only
c ---------------
      IMPLICIT    NONE
c
      include       'param-spinRL.h'
      integer*4      PP,rank,wrt
      common/ewald2/ PP
c
      integer*4      i
      real*8         intCAF, dInterpol, x
      common/intcaf/ intCAF(0:P_max,0:2*MINTPOL+1)
c                    ***************************** defined here.
c
      if(PP.ne.3) then
        write(wrt,*) ' ## P = 3 only ###'
        stop
      end if
c
      dInterpol= dfloat(MINTPOL)
c
      if(rank.eq.0) then
        write(wrt,601) PP
  601   format(/,' - interpolating the order-',i1,' charge assignment',
     *         ' function')
      end if
c
      do i= -MINTPOL, MINTPOL
      x= i/(2.d0*dInterpol)
      intCAF(0, i+MINTPOL) = 0.50d0*(0.5d0 - x)**2
      intCAF(1, i+MINTPOL) = 0.75d0 - x*x
      intCAF(2, i+MINTPOL) = 0.50d0*(0.5d0 + x)**2
      end do
c
      return
      end
c
c
c/*---------------------------------------------------------------------
      subroutine  calculate_meshift (rank,wrt)
c/*---------------------------------------------------------------------
      IMPLICIT   NONE
c
      include    'param-spinRL.h'
      integer*4   rank,wrt,i,j,k
c-----------
      real*8         meshiftx,meshifty,meshiftz,Dnx,Dny,Dnz,
     *               dmeshx,dmeshy,dmeshz
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy),
     *               meshiftz(0:meshz),Dnx(0:meshx),Dny(0:meshy),
     *               Dnz(0:meshz)
c-----------
c 
      if(rank.eq.0) then
        write(wrt,*) ' - calculating mesh-shift'
      end if
c  
      dmeshx= dfloat(meshx)
      dmeshy= dfloat(meshy)
      dmeshz= dfloat(meshz)
c
      do i= 0, meshx-1
      meshiftx(i) = i - dnint(i/dmeshx)*dmeshx
      end do
c
      do j= 0, meshy-1
      meshifty(j) = j - dnint(j/dmeshy)*dmeshy
      end do
c
      do k= 0, meshz-1
      meshiftz(k) = k - dnint(k/dmeshz)*dmeshz
      end do
c
      return
      end
c
c
c/*---------------------------------------------------------------------
      function sinc (d)
c/*---------------------------------------------------------------------
      IMPLICIT NONE
      real*8   sinc,d,epsi,c2,c4,c6,c8,pi,pid,pid2
c
      epsi =  0.1d0
      c2 = -0.1666666666667d-0
      c4 =  0.8333333333333d-2
      c6 = -0.1984126984127d-3
      c8 =  0.2755731922399d-5
c
      pi = 4.d0*DATAN(1.d0)
      pid = pi*d
c
      if (dabs(d).gt.epsi) then
         sinc = dsin(pid) / pid
      else 
         pid2 = pid*pid
         sinc = 1.d0 + pid2*( c2 + pid2*( c4 + pid2*(c6 + pid2*c8) ) )
      end if
c
      if(dabs(sinc).lt.1.d-100) sinc= 0.d0
c
      return
      end
c
c
c---------------------------------------------------------------------
       function floor (x)
c---------------------------------------------------------------------
       real*8   floor, x, xlim
c
       xlim= 100000.d0
       if(abs(x).lt.xlim) then
         floor = int(x + xlim) - xlim
       else
         write(11,*) " FLOOR: Argument too large -- Run terminated."
         call exit(1)
       end if
c
       return
       end
c
c
c----------------------------------------------------------
      subroutine magnetiz (spx,spy,spz,g,wx,wy,wz,wn,U,Uav,
     *                     wt,np,ic)
c----------------------------------------------------------
c  Magnetization
c
      implicit  none
      include  'param-spinRL.h'
      real*8    spx(np0),spy(np0),spz(np0),g,wx,wy,wz,wn,
     *          U,Uav,wt
      integer*4 np,ic,i
c
c ic< 0: reset and get instantaneous values
c ic= 0: reset
c ic= 7: get averaged value
c   one cycle must be ic= 1 -> 7 -> 0
c
      if(ic.le.0) then
        wx= 0
        wy= 0
        wz= 0
        wn= 0
        Uav= 0
        wt= 0
        if(ic.eq.0) return
      end if
c
      do 100 i= 1,np
      wx= wx +spx(i)
      wy= wy +spy(i)
      wz= wz +spz(i)
      wn= wn +1
  100 continue
c
      Uav= Uav +U      ! System total energy (not per atom)
      wt= wt +1
c
      if(ic.lt.0 .or. ic.eq.7) then
        wx= -g*wx/wn   ! in mue_b unit
        wy= -g*wy/wn
        wz= -g*wz/wn
        Uav= Uav/wt
      end if
c
      return
      end
c
c
c------------------------------------------------------------
      subroutine READ_CONF (xleng0,yleng0,zleng0,rank,suffix)
c------------------------------------------------------------
      implicit    none
      include     'param-spinRL.h'
c
      integer*4   rank,it,is,mx,my,mz,itermax,n_MCsteps,nt_P3M,
     *            ifbatch,wrt
      real*8      t,xleng,yleng,zleng,xleng0,yleng0,zleng0,
     *            rcut,pi,dt,tmax,cptot,spin2,spin3,
     *            Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,
     *            tau_diss,Temp,Tcurie,toler
      real*4      dtwr,dtwr2
      character*4 praefix*6,suffix*1,text1*40
c
      real*8        rad_fe,rad_o,elj_fe,elj_o,rcutLJ,rcutC
      common/atoms/ rad_fe,rad_o,elj_fe,elj_o,rcutLJ,rcutC
      common/parm1/ it,is
      common/parm2/ dtwr,dtwr2
      common/parm3/ t,rcut,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
      common/parm4/ mx,my,mz
      common/parm5/ n_MCsteps,nt_P3M
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,
     *              tau_diss,Temp,Tcurie
      common/itera/ toler,itermax
      common/print6/ ifbatch,wrt
c
c
      if(rank.eq.0) then
        write(wrt,*) ' read: ',praefixs//'_config.START'//suffix
        write(wrt,*) ' write: ',praefixc//'.13'//suffix
        write(wrt,*) 'READ_CONF: Parameter read... start'
      end if
c
      OPEN (unit=08,file=praefixs//'_config.START'//suffix,
     *                                   status='old',form='formatted')
c
c* ------------------------------------------------------
      read (8,'(a40,a6)') text1, praefix   ! String of identification
      read (8,'(a40,f20.0)') text1,cptot   ! Maximum wall time for each run 
      read (8,'(a40,f20.0)') text1,tmax    ! Maximum physical time (ps)
      read (8,'(a40,f20.0)') text1,dt      ! Time step (ps)
      read (8,'(a40,f20.0)') text1,dtwr    ! Energy write-out interval
      read (8,'(a40,f20.0)') text1,dtwr2   ! Data output interval 
      read (8,'(a40,f20.0)') text1,xleng0  ! Length x of unit cell (Ang)
      read (8,'(a40,f20.0)') text1,yleng0  ! Length y of unit cell 
      read (8,'(a40,f20.0)') text1,zleng0  ! Length z of unit cell 
      read (8,'(a40,i10)')   text1,mx      ! Number of domains
      read (8,'(a40,i10)')   text1,my      !
      read (8,'(a40,i10)')   text1,mz      !
      read (8,'(a40,i10)')   text1,n_MCsteps ! Maximum iteration count of MC
      read (8,'(a40,f20.0)') text1,rcut    ! Cutoff radius of spin interaction
      read (8,'(a40,f20.0)') text1,rcutC   ! Cutoff radius of ES interaction
      read (8,'(a40,i10)')   text1,itermax ! Maximum iteration of nonlinear solver 
      read (8,'(a40,f20.0)') text1,toler   ! convergence or iteration tolerance
      read (8,'(a40,f20.0)') text1,Jaa     ! Exchange integral A-A (10meV)
      read (8,'(a40,f20.0)') text1,Jbb     ! Exchange integral B-B 
      read (8,'(a40,f20.0)') text1,Jab     ! Exchange integral A-B 
      read (8,'(a40,f20.0)') text1,Bapx    ! Magnetic field Bx (100gauss)
      read (8,'(a40,f20.0)') text1,Bapy    ! Magnetic field By
      read (8,'(a40,f20.0)') text1,Bapz    ! Magnetic field Bz
      read (8,'(a40,f20.0)') text1,tau_b   ! Period of microwave (ps)
      read (8,'(a40,f20.0)') text1,tau_diss! Relaxation time (ps)
      read (8,'(a40,f20.0)') text1,Temp    !  Temperature in Kelvin
      read (8,'(a40,f20.0)') text1,Tcurie  !  Curie Temperature in Kelvin
      read (8,'(a40,f20.0)') text1,rcutLJ  !  rcut of LJ force
c
      close (8)
c
      return
      end
c
c
c-----------------------------------------------------------------------
      subroutine Temp_fix (spx,spy,spz,sp2,aTsz,Tsz0,spec,site,np,
     *                     ifcmp,nstep_mcT)
c-----------------------------------------------------------------------
      implicit none
c
      include    'mpif.h'
      include    'param-spinRL.h'
c
      real*8     spx(np0),spy(np0),spz(np0),sp2(np0),
     *           aTsz,Tsz0,ss,sp,sc,alpha,pi,th,ph,ranff,th0,
     *           sA3,sB3,sB2
      integer*4  spec(np0),site(np0),np,ifcmp,nstep_mcT,i_mc,i,j,
     *           nA3,nB3,nB2
c
c ---------------------------------
      if(ifcmp.eq.1) go to 2000
c ---------------------------------
c
c* Radomization...
c
      i_mc= 0
      pi= 4*datan(1.d0)
c
 1000 i_mc= i_mc +1
      if(i_mc.gt.nstep_mcT) return   ! No operation if nstep_mcT= 0
c
      j= mod(int(np*ranff(0.)+0.00001),np) +1
      th= pi*ranff(0.)
      ph= 2*pi*ranff(0.)
c
      th0= dacos(spz(j)/sp2(j))
      th = 0.7d0*th0 +0.3d0*th    ! weighted average 
c 
cc    spx(j)= sp2(j)*dsin(th)*dcos(ph)
cc    spy(j)= sp2(j)*dsin(th)*dsin(ph)
      spz(j)= sp2(j)*dcos(th)
c
      ss= sqrt(sp2(j)**2 -spz(j)**2)
      sp= sqrt(spx(j)**2 +spy(j)**2)
c
      spx(j)= ss *spx(j)/sp       ! keep the rotation phase
      spy(j)= ss *spy(j)/sp
c
      go to 1000
c ----------------------------
c
 2000 continue
c
c* Compression toward its axis ...
c
      sA3= 0     ! A site,3+
      sB3= 0     ! B site,3+
      sB2= 0     ! B site,2+
      nA3= 0
      nB3= 0
      nB2= 0
c
      do 300 i= 1,np
      if(site(i).eq.1) then
        sA3= sA3 +spz(i)
        nA3= nA3 +1
      else
        if(spec(i).eq.1) then
          sB3= sB3 +spz(i)
          nB3= nB3 +1
        else
          sB2= sB2 +spz(i)
          nB2= nB2 +1
        end if
      end if
  300 continue
c
      sA3= sA3/nA3   ! axis of distribution
      sB3= sB3/nB3
      sB2= sB2/nB2
c
c* compress the length  -S.......0.......S
c                        I     s         I
c               -S+alp*(s+S)         S-alp*(S-s)
c               =(-S)*(1-alp)+alp*s  =S*(1-alp)+alp*s
c
      alpha= Tsz0/aTsz
c
      do 400 i= 1,np
      if(site(i).eq.1) then
        if(sA3.lt.0.) sc= -sp2(i)
        if(sA3.gt.0.) sc=  sp2(i)
      else
        if(spec(i).eq.1) then
          if(sB3.lt.0.) sc= -sp2(i)
          if(sB3.gt.0.) sc=  sp2(i)
        else
          if(sB2.lt.0.) sc= -sp2(i)
          if(sB2.gt.0.) sc=  sp2(i)
        end if
      end if
c
      spz(i)= sc*(1.d0 -alpha) + alpha*spz(i)
c
      ss= sqrt(sp2(i)**2 -spz(i)**2)
      sp= sqrt(spx(i)**2 +spy(i)**2)
c
      spx(i)= ss *spx(i)/sp
      spy(i)= ss *spy(i)/sp
  400 continue
c
      return
      end
c
c
c-------------------------------------------------------------
      subroutine plot_spin (x,y,z,spx,spy,spz,spec,site,np)
c-------------------------------------------------------------
      include  'param-spinRL.h'
c
      real*8     x(np0),y(np0),z(np0),spx(np0),spy(np0),spz(np0)
      integer*4  spec(np0),site(np0)
      real*8     t,xleng,yleng,zleng,rcut,pi,dt,tmax,cptot,
     *           spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,tau_diss,
     *           Temp,Tcurie
      common/parm3/ t,rcut,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,
     *              tau_diss,Temp,Tcurie
c
      CHARACTER*8    LABEL,cax*1
      COMMON/HEADR1/ LABEL,cdate
      COMMON/HEADR2/ time
c
      HH= 0.7
      CALL SYMBOL ( 0.5,17.0,HH,LABEL,0.,8)
      CALL SYMBOL ( 5.0,17.0,HH,praefixc,0.,24)
      CALL SYMBOL (15.0,17.0,HH,'Spin dynamics (MD) ',0.,20)
c
      CALL SYMBOL ( 0.0,1.0,HH,cdate, 0.,9)
      CALL SYMBOL (19.0,1.0,HH,'Time=', 0.,5)
      CALL VALUES (20.5,1.0,HH,time,0.,101)
c
      spin24= spin2
      spin34= spin3
      Tcurie4= Tcurie
      CALL SYMBOL (0.5,16.0,HH,'Spin(A)=', 0.,8)
      CALL VALUES (3.5,16.0,HH,spin24,0.,101)
      CALL SYMBOL (0.5,15.2,HH,'Spin(B)=', 0.,8)
      CALL VALUES (3.5,15.2,HH,spin34,0.,101)
      CALL SYMBOL (0.5,14.4,HH,'Tcurie=', 0.,7)
      CALL VALUES (3.5,14.4,HH,Tcurie4,0.,101)
c
      qaa4= 10 * Jaa  ! in 10 meV
      qbb4= 10 * Jbb
      qab4= 10 * Jab
      CALL SYMBOL (7.5,16.0,HH,'Jaa=', 0.,4)
      CALL VALUES (9.0,16.0,HH,qaa4,0.,101)
      CALL SYMBOL (7.5,15.2,HH,'Jbb=', 0.,4)
      CALL VALUES (9.0,15.2,HH,qbb4,0.,101)
      CALL SYMBOL (7.5,14.4,HH,'Jab=', 0.,4)
      CALL VALUES (9.0,14.4,HH,qab4,0.,101)
c
      Temp4= Temp
      Bext4= sqrt(Bapx**2 +Bapy**2 +Bapz**2)
      CALL SYMBOL (13.0,16.0,HH,'Temp=', 0.,5)
      CALL VALUES (15.0,16.0,HH,Temp4,0.,101)
      CALL SYMBOL (13.0,15.2,HH,'Bapp=', 0.,5)
      CALL VALUES (15.0,15.2,HH,Bext4,0.,101)
c
      xleng4= xleng
      dt4 = dt
      tau4= tau_b
      tau_diss4= tau_diss
c     CALL SYMBOL (19.0,16.0,HH,'Leng=', 0.,5)
c     CALL VALUES (21.0,16.0,HH,xleng4,0.,101)
      CALL SYMBOL (19.0,16.0,HH,'tau_di=', 0.,7)
      CALL VALUES (21.0,16.0,HH,tau_diss4,0.,101)
      CALL SYMBOL (19.0,15.2,HH,'np=', 0.,3)
      CALL VALUES (21.0,15.2,HH,float(np),0.,101)
      CALL SYMBOL (19.0,14.4,HH,'dt=', 0.,3)
      CALL VALUES (21.0,14.4,HH,dt4,0.,101)
      CALL SYMBOL (19.0,13.6,HH,'tau_b=', 0.,6)
      CALL VALUES (21.0,13.6,HH,tau4,0.,101)
c
      FSIZE= 4.
      HL=  11.
      VD=   6.
      phi= -60.
      tht=  15.
c
      pi=  4.*atan(1.0)
      pha= pi*phi/180.
      tha= pi*tht/180.
c
      cph= cos(pha)
      sph= sin(pha)
      cth= cos(tha)
      sth= sin(tha)
c
      xmax= 0.5*xleng
      ymax= 0.5*yleng
      zmax= 0.5*zleng
c
      xp= xmax*cph -ymax*sph
      yp= xmax*sph +ymax*cph
      zp= zmax
c
      ypp=  yp
      zpp= -xp*sth +zp*cth
c
      rmax1= sqrt(ypp**2 +zpp**2)
      ps= fsize/rmax1
c
c**********************
C*  Draw Axes.        *
c**********************
c
      do 100 i= 1,3
      if(i.eq.1) then
         x1= xmax
         y1= 0.
         z1= 0.
         cax='X'
      else if(i.eq.2) then
         x1= 0.
         y1= ymax
         z1= 0.
         cax='Y'
      else if(i.eq.3) then
         x1= 0.
         y1= 0.
         z1= zmax
         cax='Z'
      end if
c
      xp= x1*cph -y1*sph
      yp= x1*sph +y1*cph
      zp= z1
c
      ypp= yp
      zpp= -xp*sth +zp*cth
c
      xx= ps*ypp  +HL
      yy= ps*zpp  +VD
      CALL PLOT (HL,VD,3)
      CALL PLOT (xx,yy,2)
c
      CALL SYMBOL (xx-0.7,yy-0.5,HH,cax,0.,1)
  100 continue
c
c-------------------
      dd= 0.2
      ff= 0.25
c-------------------
      do 400 i= 1,np
      xp= (x(i) -ff*spx(i))*cph -(y(i) -ff*spy(i))*sph
      yp= (x(i) -ff*spx(i))*sph +(y(i) -ff*spy(i))*cph
      zp=  z(i) -ff*spz(i) 
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
      xx1= ps*ypp +HL
      yy1= ps*zpp +VD
      if(xx1.lt.0. .or. xx1.gt.20.) go to 400
      if(yy1.lt.0. .or. yy1.gt.18.) go to 400
c
      xp= (x(i) +ff*spx(i))*cph -(y(i) +ff*spy(i))*sph
      yp= (x(i) +ff*spx(i))*sph +(y(i) +ff*spy(i))*cph
      zp=  z(i) +ff*spz(i) 
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
      xx2= ps*ypp +HL
      yy2= ps*zpp +VD
c
      if(spec(i).eq.1) then      !* Fe(3+)
        if(site(i).eq.1) then    ! A-site
          call newcolor (3,0.,1.,0.)  ! 3-color (r,g,b)
          call circle (xx1-0.01,yy1-0.05,dd,2)
        else                     ! B-site
          call newcolor (3,0.,0.,1.)  ! 3-color (r,g,b)
          call circle (xx1-0.01,yy1-0.05,dd,2)
        end if
      end if
      if(spec(i).eq.2) then     !* Fe(2+)
        call newcolor (3,1.,0.,0.)  ! 3-color (r,g,b)
        call triang (xx1-0.01,yy1-0.05,dd,2)
      end if
c
      call plot (xx1,yy1,3)
      call plot (xx2,yy2,2)
      call plot (xx2,yy2,3)
  400 continue
c
      call newcolor (0,1.,0.,0.)  ! gray scale
C---------------------
      CALL CHART
C---------------------
C
      return
      end
c
c
c-----------------------------------------------------------
      subroutine distr_spin (spx,spy,spz,spec,site,np)
c-----------------------------------------------------------
      include   'param-spinRL.h'
c
      real*8     spx(np0),spy(np0),spz(np0),ss
      real*4     ang1x(101),ang2x(101),ang1y(101),ang2y(101),
     *           ang1z(101),ang2z(101),hha(101),
     *           angax(101),angay(101),angaz(101)
      integer*4  spec(np0),site(np0),np
c
      do 100 k= 1,101
      ang1x(k)= 0.
      ang2x(k)= 0.
      ang1y(k)= 0.
      ang2y(k)= 0.
      ang1z(k)= 0.
      ang2z(k)= 0.
      angax(k)= 0.
      angay(k)= 0.
      angaz(k)= 0.
      hha(k)= (k-51)/50.
  100 continue
c
      do 200 i= 1,np
      ss= sqrt(spx(i)**2+spy(i)**2+spz(i)**2)
c
      costhx= spx(i)/ss
      costhy= spy(i)/ss
      costhz= spz(i)/ss
c
      ir= 50.*(costhx +1.) +1.01
      if(ir.gt.0 .and. ir.lt.101) then
        if(site(i).eq.1) then
          angax(ir)= angax(ir) +1.                    ! A:3+
        else
          if(spec(i).eq.1) ang1x(ir)= ang1x(ir) +1.   ! B:3+
          if(spec(i).eq.2) ang2x(ir)= ang2x(ir) +1.   ! B:2+
        end if
      end if
c
      ir= 50.*(costhy +1.) +1.01
      if(ir.gt.0 .and. ir.lt.101) then
        if(site(i).eq.1) then
          angay(ir)= angay(ir) +1.
        else
          if(spec(i).eq.1) ang1y(ir)= ang1y(ir) +1.
          if(spec(i).eq.2) ang2y(ir)= ang2y(ir) +1.
        end if
      end if
c
      ir= 50.*(costhz +1.) +1.01
      if(ir.gt.0 .and. ir.lt.101) then
        if(site(i).eq.1) then
          angaz(ir)= angaz(ir) +1.
        else
          if(spec(i).eq.1) ang1z(ir)= ang1z(ir) +1.
          if(spec(i).eq.2) ang2z(ir)= ang2z(ir) +1.
        end if
      end if
  200 continue
c
      ILN= -1
      ILG= 2
      call lplmax (ang1x,famax1,famin,101)
      call lplmax (ang2x,famax2,famin,101)
      famax= amax1(famax1,famax2)
      nxtick= 3
      nytick= 3
      call hplot1 (2,2,101,hha,ang1x,famax,famin,ILN,nxtick,nytick,
     *               '        ',8,'cos(thx)',8,' (3 B   ',8,0)
      call hplot1 (2,2,101,hha,angax,famax,famin,ILN,nxtick,nytick,
     *               '        ',8,'cos(thx)',8,' (3A    ',8,1)
      call hplot1 (2,3,101,hha,ang2x,famax,famin,ILN,nxtick,nytick,
     *               '        ',8,'cos(thx)',8,' (2B)   ',8,0)
c
      call lplmax (ang1y,famax1,famin,101)
      call lplmax (ang2y,famax2,famin,101)
      famax= amax1(famax1,famax2)
      call hplot1 (3,2,101,hha,ang1y,famax,famin,ILN,nxtick,nytick,
     *               '        ',8,'cos(thy)',8,' (3 B   ',8,0)
      call hplot1 (3,2,101,hha,angay,famax,famin,ILN,nxtick,nytick,
     *               '        ',8,'cos(thy)',8,' (3A    ',8,1)
      call hplot1 (3,3,101,hha,ang2y,famax,famin,ILN,nxtick,nytick,
     *               '        ',8,'cos(thy)',8,' (2B)   ',8,0)
C---------------------
      CALL CHART
C---------------------
c
      call lplmax (ang1z,famax1,famin,101)
      call lplmax (ang2z,famax2,famin,101)
      famax= amax1(famax1,famax2)
      call hplot1 (2,2,101,hha,ang1z,famax,famin,ILN,nxtick,nytick,
     *               '        ',8,'cos(thz)',8,' (3 B   ',8,0)
      call hplot1 (2,2,101,hha,angaz,famax,famin,ILN,nxtick,nytick,
     *               '        ',8,'cos(thz)',8,' (3A    ',8,1)
      call hplot1 (2,3,101,hha,ang2z,famax,famin,ILN,nxtick,nytick,
     *               '        ',8,'cos(thz)',8,' (2B)   ',8,0)
C---------------------
      CALL CHART
C---------------------
      return
      end
c
c
c------------------------------------------------------
      subroutine plot_disp (axis,freq,modes,plot_ch)
c------------------------------------------------------
      real*4   axis(100),freq(100)
      character*8  plot_ch
      integer*4  ifbatch,wrt
      common/print6/ ifbatch,wrt
c
      ILN= 1
      nxtick= 4
      nytick= 4
      call lplmax (freq,emax,emin,modes)
        write(wrt,*) 'plot_disp: max,min=',emax,emin
c
      call lplot1 (1,1,modes,axis,freq,emax,emin,ILN,nxtick,
     *             nytick,plot_ch,8,'wave num',8,'        ',8,0)
c---------------------
      CALL CHART
c---------------------
      return
      end
c
c
c------------------------------------------------------
      subroutine lplots (suffix)
c------------------------------------------------------
      implicit      none
c
      include       'param-spinRL.h'
      real*4        spinx,spinz,spin7,bextx,bextz,magx,magy,magz,
     *              Usys,conv,aitr,psdt,Tfix,Uss,Usb,Tsx,Tsy,Tsz,
     *              sum_mb,U_fe,U_o,ds_fe,ds_o,fdt4,vdt4,idt4,timeh,
     *              cosd,sind,chi_real,chi_imag,B00,Bmw,
     *              mue_b,hh,av_mz(nhs),amz,ss,s2,
     *              emax,emin,emax1,emax2,emax3,emin1,emin2,emin3
      common/ehist/ spinx(nhs),spinz(nhs),spin7(nhs),
     *              bextx(nhs),bextz(nhs),magx(nhs),magy(nhs),magz(nhs),
     *              Usys(nhs),conv(nhs),aitr(nhs),psdt(nhs),Tfix(nhs),
     *              Uss(nhs),Usb(nhs),Tsx(nhs,3),Tsy(nhs,3),Tsz(nhs,3),
     *              sum_mb(nhs),U_fe(nhs),U_o(nhs),ds_fe(nhs),ds_o(nhs),
     *              fdt4(nhs),vdt4(nhs),idt4(nhs),timeh(nhs)
c
      real*8        spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,
     *              tau_diss,Temp,Tcurie
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,
     *              tau_diss,Temp,Tcurie
      real*8        t,xleng,yleng,zleng,rcut,pi,dt,tmax,cptot
      common/parm3/ t,rcut,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
c
      CHARACTER*8    LABEL,suffix*1
      integer*4      i,k,it,is,is0,ns,
     *               ILN,ILG,nxtick,nytick
      COMMON/HEADR1/ LABEL,CDATE
      common/parm1/  it,is
      integer*4      ifbatch,wrt
      common/print6/ ifbatch,wrt
c
      ILN= 1
      ILG= 2
      hh= 0.7
      CALL SYMBOL (1.0,18.5,hh,'Spin dynamics (MD)  ',0.,20)
c
c-------------------------------------------------
c* One-period averaged Mz(t)= Mz - <Mz>
c   40 data points are generated per period
c-------------------------------------------------
c
      is0= av_start*is  ! start averaging at this step
c
      do 103 k= 1,is
      av_mz(k)= magz(k)
  103 continue
c
      do 100 k= 40,is-40
      ss= 0
      ns= 0
c
      do 130 i= k-39,k+40
      ss= ss +magz(i)
      ns= ns +1
  130 continue
c
      av_mz(k)= ss/ns
  100 continue
c
c
      ss= 0
      s2= 0
      ns= 0
c
      do 200 k= is0,is-40
      ss= ss +(magz(k) -av_mz(k))*bextz(k)   ! magz= -g*<sz>
      s2= s2 +(magz(k) -av_mz(k))**2         ! <dm**2>
      ns= ns +1
  200 continue
c
      mue_b= 9.27410e-21   ! e*hbar/2mc
      B00= 100             ! gauss
      Bmw= B00*sqrt(Bapx**2 +Bapy**2 +Bapz**2)
c
      amz= sqrt(2*s2/ns)
c
      if(Bmw.gt.1.d-5) then
        chi_real= mue_b*amz/Bmw
        chi_imag= 2*(mue_b*B00*ss/ns)/Bmw**2
c
        write(wrt,731) is,amz,chi_real,chi_imag
  731   format('Magnetic susceptibility...',
     *         '   is, <Mz> =',i5,f8.3,/,
     *         '   chi_r, chi_i=',1p2e11.3,/)
      else
        write(wrt,*) 'Magnetic susceptibility is not'
        write(wrt,*) '  calculated as Bmw = 0 ......'
      end if
c
c
      call lplmax (magx,emax1,emin1,is)
      call lplmax (magy,emax2,emin2,is)
      call lplmax (magz,emax3,emin3,is)
      emax = amax1(emax1,emax2,emax3,-emin1,-emin2,-emin3)
      emin = -emax
      nxtick= 3
      nytick= 3
      call lplot1 (2,4,is,timeh,magx,emax,emin,ILN,nxtick,nytick,
     *             'Magnet-x',8,'        ',8,'        ',8,0)
      call lplot1 (2,5,is,timeh,magy,emax,emin,ILN,nxtick,nytick,
     *             'Magnet-y',8,'        ',8,'        ',8,0)
      call lplot1 (2,6,is,timeh,magz,emax3,emin3,ILN,nxtick,nytick,
     *             'Magnet-z',8,' time   ',8,'        ',8,0)
c
c     call lplmax (spinx,emax1,emin1,is)
      call lplmax (spinz,emax2,emin2,is)
      emax = amax1(emax1,emax2,-emin1,-emin2)
      emin = -emax
c     call lplot1 (3,4,is,timeh,spinx,emax,emin,ILN,nxtick,nytick,
c    *             'Spin-x.7',8,'        ',8,'        ',8,0)
      call lplot1 (3,5,is,timeh,spinz,emax,emin,ILN,nxtick,nytick,
     *             'Spin-z.7',8,'        ',8,'        ',8,0)
c
      call lplmax (spin7,emax,emin,is)
      emax= 1.2*emax
      call lplot1 (3,6,is,timeh,spin7,emax,0.,ILN,nxtick,nytick,
     *             'Spin-7  ',8,' time   ',8,'        ',8,0)
c------------------------
      CALL CHART
c------------------------
c
      call lplmax (Usys,emax,emin,is)
      call lplot1 (2,4,is,timeh,Usys,emax,emin,ILN,nxtick,nytick,
     *             'Usys    ',8,'        ',8,'        ',8,0)
c
      call lplmax (bextx,emax,emin,is)
      emax = amax1(emax,-emin)
      emin = -emax
      call lplot1 (2,5,is,timeh,bextx,emax,emin,ILN,nxtick,nytick,
     *             'bextx.7 ',8,'        ',8,'        ',8,0)
c
      call lplmax (bextz,emax,emin,is)
      emax = amax1(emax,-emin)
      emin = -emax
      call lplot1 (2,6,is,timeh,bextz,emax,emin,ILN,nxtick,nytick,
     *             'bextz.7 ',8,'  time  ',8,'        ',8,0)
c
      call lplmax (Uss,emax,emin,is)
      call lplot1 (3,4,is,timeh,Uss,emax,emin,ILN,nxtick,nytick,
     *             ' Us*s   ',8,'        ',8,'        ',8,0)
c
      call lplmax (Usb,emax,emin,is)
      call lplot1 (3,5,is,timeh,Usb,emax,emin,ILN,nxtick,nytick,
     *             ' Us*B   ',8,'        ',8,'        ',8,0)
c
      CALL SYMBOL (18.5,3.8,hh,'del.Mz= ', 0.,8)
      CALL VALUES (22.0,3.8,hh,amz,0.,101)
c
      CALL SYMBOL (18.5,2.8,hh,'Re_chi =', 0.,8)
      CALL VALUES (22.0,2.8,hh,chi_real,0.,101)
      CALL SYMBOL (18.5,2.1,hh,'Im_chi =', 0.,8)
      CALL VALUES (22.0,2.1,hh,chi_imag,0.,101)
c------------------------
      CALL CHART
c------------------------
c
      call lplmax (ds_fe,emax,emin,is)
      emax = amax1(emax,-emin)
      emin = -emax
      call lplot1 (2,4,is,timeh,ds_fe,emax,emin,ILN,nxtick,nytick,
     *             ' ds_Fe  ',8,'        ',8,'        ',8,0)
c
      call lplmax (ds_o,emax,emin,is)
      emax = amax1(emax,-emin)
      emin = -emax
      call lplot1 (2,5,is,timeh,ds_o,emax,emin,ILN,nxtick,nytick,
     *             ' ds_O   ',8,'        ',8,'        ',8,0)
c
      call lplmax (U_o,emax,emin,is)
      call lplot1 (3,4,is,timeh,U_o,emax,emin,ILN,nxtick,nytick,
     *             ' Kin_O  ',8,'        ',8,'        ',8,0)
c
      call lplmax (U_fe,emax,emin,is)
      call lplot1 (3,5,is,timeh,U_fe,emax,emin,ILN,nxtick,nytick,
     *             ' Kin_Fe ',8,'        ',8,'        ',8,0)
c
      call lplmax (fdt4,emax,emin,is)
      call lplot1 (2,6,is,timeh,fdt4,emax,emin,ILN,nxtick,nytick,
     *             'F*dt/mvo',8,'  time  ',8,'        ',8,0)
c
      call lplmax (vdt4,emax,emin,is)
      call lplot1 (3,6,is,timeh,vdt4,emax,emin,ILN,nxtick,nytick,
     *             ' V*dt   ',8,'  time  ',8,'        ',8,0)
c------------------------
      CALL CHART
c------------------------
c
      return
      end
c
c
c------------------------------------------------------
      subroutine rehist (rank,suffix)
c------------------------------------------------------
      implicit      none
      include       'param-spinRL.h'
c
      integer*4     rank,i,j,it,is,iss,iwa,iwb,ifbatch,wrt
      common/parm1/ it,is
      common/imemo/ iwa,iwb
c
      real*4        array,dtwr,dtwr2
      common/ehist/ array(nhs,33)
      common/parm2/ dtwr,dtwr2
      real*8        t,xleng,yleng,zleng,rcut,pi,dt,tmax,cptot
      common/parm3/ t,rcut,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
      common/print6/ ifbatch,wrt
      character     suffix*1
c
      do 100 j= 1,33
      iss= 0
c
      do 100 i= 1,is,2
      iss= iss +1
      array(iss,j)= 0.5*(array(i,j) +array(i+1,j))
  100 continue
c
      is= iss
      dtwr= 2*dtwr
      iwa= 0 ! t/dtwr
c
      if(rank.eq.0) then
        write(wrt,*) ' Rehist is called: it, is=',it,is
      end if
c
      return
      end
c
c
c---------------------------------------------------------------
      subroutine ggauss
c---------------------------------------------------------------
      common/gaus1/ fv(51),vv0,dv
c
      FV(1)=0.0
C
      VV0= -3.
      DV= 2.*ABS(VV0)/50.0
C
      VV= VV0
      DO 100 J=1,50
      S=0.0
      NS=1000
      K2=NS/2
      SDV=DV/FLOAT(NS)
C
      DO 130 K=1,K2
      VV=VV +2.0*SDV
      S=S +4.0*FUN(VV-SDV) +2.0*FUN(VV)
  130 CONTINUE
      S= (S +4.0*FUN(VV+SDV) +FUN(VV+2.0*SDV))*SDV/3.0
      FV(J+1)= FV(J)+S
  100 CONTINUE
C
      DO 200 I=1,51
      FV(I)=FV(I)/FV(51)
  200 CONTINUE
C
      return
      end
c
c
c---------------------------------------------------------------
      real*8 function dgaus2 (vmax)
c---------------------------------------------------------------
      real*8   eps,ranff,vmax
      common/gaus1/ fv(51),vv0,dv
c
      eps= ranff(0.)
      DO 100 K=1,51
      K2=K
      IF(FV(K).GT.eps) GO TO 200
  100 CONTINUE
C
  200 Y1= FV(K2-1)
      Y2= FV(K2)
      X2= (EPS-Y2)/(Y2-Y1)+K2
      dgaus2= vmax*(VV0 +DV*(X2-1.0))
c
      return
      end
c
c
c---------------------------------------------------------------
      function fun (v)
c---------------------------------------------------------------
      fun= exp(-v**2/2.)
c
      return
      end
c
c
c------------------------------------------------------
      subroutine averg1 (q,qav,is)
c------------------------------------------------------
      dimension  q(5000)
c
      qav= 0.
c
      do 100 i= is-9,is
      qav= qav +q(i)
  100 continue
c
      qav= qav/10.
c
      return
      end
c
c
c------------------------------------------------------
      subroutine clocks (cputime,walltime)
c------------------------------------------------------
c* Get the cpu and elapsed times on each node.
c
      include 'mpif.h'
c
      real*4  cputime,ct,tm(2)
      real*8  walltime,walltime0
      logical first_clk
      common/saveif2/ walltime0,first_clk
c               - defined as .true. in /block data/
c
      if(first_clk) then
         walltime0= mpi_wtime()
         first_clk= .false.
      end if
c
c     ct = etime(tm)   ! Pentium and Unix general.
      cputime = 0      ! tm(1)
      walltime= mpi_wtime() -walltime0  ! sec
c
      return
      end
c
c
c------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
c------------------------------------------------------
      include    'param-spinRL.h'
      real*4     f(nhs),fmax,fmin
      integer*4  is
c
      fmax= -1.e10
      fmin=  1.e10
c
      do 100 i= 1,is
      fmax= amax1(fmax,f(i))
      fmin= amin1(fmin,f(i))
  100 continue
c
      return
      end
c
c
c-------------------------------------------------
       subroutine circle (x,y,d,ic)
c-------------------------------------------------
c*  Open circle centered at (x,y) /or outer edge.
c
      write(77,*) " 3.0 setlinewidth"
c
      pi= 3.1415927
      nc= 13
      dth= 2.*pi/nc
      a= d/2.
c
      x0= x +a
      y0= y
      call plot (x0,y0,3)
c
      do 100 j= 1,nc
      th= dth*j
c
      x1= x +a*cos(th)
      y1= y +a*sin(th)
c
      call plot (x1,y1,2)
  100 continue
c
      call plot (x1,y1,3)
      write(77,*) " 1.0 setlinewidth"
c
      if(ic.eq.1) return
c------------------------------------
c*  Filled circle centered at (x,y).
c------------------------------------
c
      write(77,*) " 3.0 setlinewidth"
c
      nc= 5
      dth= pi/(2*nc +1)
c
      do 300 j= -nc,nc
      th= 0.5*pi +dth*j
c
      x1= x +a*cos(th)
      y1= y +a*sin(th)
c
      x2= x1
      y2= 2.*y -y1
c
      call plot (x1,y1,3)
      call plot (x2,y2,2)
  300 continue
c
      call plot (x2,y2,3)
      write(77,*) " 1.0 setlinewidth"
c
      return
      end
c
c
c-------------------------------------------------
        subroutine triang (x,y,d,ic)
c-------------------------------------------------
c*  Open triangle centered at (x,y) /or outer edge.
c
      write(77,*) " 3.0 setlinewidth"
c
      a= d/2.
      b= a/1.732
      c= 2.*b
c
      call plot (x-a,y-b,3)
      call plot (  x,y+c,2)
      call plot (x+a,y-b,2)
      call plot (x-a,y-b,2)
      call plot (x+a,y-b,3)
c
      write(77,*) " 1.0 setlinewidth"
      if(ic.eq.1) return
c
c------------------------
c*  Fill the triangle.
c------------------------
      nc=7
      nch= (nc+1)/2
      dx= d/nc
      y2= y -b
c
      do 100 j= 1,nc
      x1= x-a +dx*(j-0.5)
      if(j.le.nch) then
         y1= y +1.732*a*(j-0.5)/nch
      else
         y1= y +1.732*a*(nc-j)/nch 
      end if
c
      call plot (x1,y1,3)
      call plot (x1,y2,2)
  100 continue
c
      call plot (x1,y2,3)
      write(77,*) " 1.0 setlinewidth"
c
      return
      end
c
c
c------------------------------------------------
      integer*4 function iwrta (t,dtwr)
c------------------------------------------------
      real*8  t
      common/imemo/ iwa,iwb
c
      iw= t/dtwr
      if(iw.gt.iwa) then
        iwa= iw
        iwrta= 0
      else
        iwrta= 1
      end if
c
      return
      end
c
c
c------------------------------------------------
      integer*4 function iwrtb (t,dtwr)
c------------------------------------------------
      real*8  t
      common/imemo/ iwa,iwb
c
      iw= t/dtwr 
      if(iw.gt.iwb) then
        iwb= iw
        iwrtb= 0
      else
        iwrtb= 1
      end if
c
      return
      end
c
c
c------------------------------------------------
      block data
c------------------------------------------------
      real*8  walltime0
      logical first_shf,first_clk
c
      common/saveif2/ walltime0,first_shf,first_clk
      common/ranfff/ ir,ir0
c
      data  first_shf/.true./, first_clk/.true./
      data  ir/17331/
      end
c
c
c------------------------------------------------
      real*8 function ranff (x)
c------------------------------------------------
c*  ranf= (0,1)
c
      common/ranfff/ ir,ir0
C
      REAL*8     INVM
      PARAMETER  (MASK=2**30+(2**30-1),INVM= 0.5D0**31)
      PARAMETER  (LAMBDA=48828125)
C
      IR= IAND( LAMBDA*IR, MASK)
      ranff= IR*INVM
c
c     ask= 371597.
c     ambda= sqrt(ask)
c     qq= 0.3713*ask
c
c     ir= amod( ambda*ir +qq, ask)
c     ranff= ir/ask
c
      return
      end
c
c
C-----------------------------------------------------------------------
      SUBROUTINE LPLOT1 (IX,IY,NPT1,X,Y,YMAX,YMIN,IL,nxtick,
     *                   nytick,LAB1,N1,LAB2,N2,LAB3,N3,iskip)
C-----------------------------------------------------------------------
C  <<Warning>>  Order and number of arguments /LPLOT/ have been changed.
C               Also, X (time) is defined for all range.
C               Date: 5/18/96 at MIT.
C***********************************************************************
C   IL=1................ LINEAR PLOT OF (X,Y)
C   IL=2................ LOG10 PLOT OF (X,LOG Y)
C***********************************************************************
C
      DIMENSION  X(5000),Y(5000),U(5000),V(5000)
      DIMENSION  XCM(6),YCM(7),PL(6),PR(6),QL(7),QR(7)
      real*4     ymax,ymin
C
      CHARACTER*8    LAB1,LAB2,LAB3,LABEL
      COMMON/HEADR1/ LABEL,cdate
      COMMON/HEADR2/ time,xleng
      COMMON/PPLCOM/ NFINE,PL1(10),PR1(10),QL1(10),QR1(10),
     *               XMIN1(10),XMAX1(10),YMIN1(10),YMAX1(10)
C
      DATA  XCM/21.0, 2*10.00, 3*6.00/,
     *      YCM/15.0, 2*6.80, 4*4.00/,
     *      PL/2.0,  2.0,14.0, 2.0,9.0,16.0/,
     *      QL/2.3, 10.5,2.3, 14.0,9.5,5.0,0.5/
C
      call set_width (1.0)
c
      IPLOT=1
      GO TO 1
C
C-----------------------------------------------------------------------
      ENTRY HPLOT1 (IX,IY,NPT1,X,Y,ymax,ymin,IL,LAB1,N1,LAB2,N2,LAB3,N3)
C-----------------------------------------------------------------------
      IPLOT=2
C
    1 NPT= NPT1
      ISC= 1
C
      DO 5 I=1,6
    5 PR(I)= PL(I) +XCM(I)
C
      DO 6 J=1,7
    6 QR(J)= QL(J) +YCM(J)
C
C                 ******************************************************
C*                **  MAKE A COPY BEFORE THE TOP-LEFT FRAME IS DRAWN. **
C                 ******************************************************
      HH= 0.70
      I1= IABS(IX)
      J1= IABS(IY)
      IF(I1.GE.3) GO TO 10
      IF(J1.EQ.3.OR.J1.GE.5) GO TO 10
C                                              ************************
C                                              ** LABEL OF THE PAGE. **
C                                              ************************
      CALL SYMBOL (20.3,0.1,HH,'T=',0.,2)
      CALL VALUES (21.5,0.1,HH,time,0.,101)
C
   10 CONTINUE
C
      DO 23 I=1,NPT
   23 U(I)= X(I)
      XMAX= U(NPT)
      XMIN= U(1)
C                             ************************************
C                             ** THREE-POINT AVERAGE IF IL > 0  **
C                             ************************************
      IF(IL.GT.0) THEN
        V(1)=   Y(1)
        V(NPT)= Y(NPT)
        v(npt-1)= y(npt-1)
c
        DO 37 I=2,NPT-2
   37   V(I)= 0.33333*(Y(I-1)+Y(I)+Y(I+1))
      ELSE
        DO 38 I=1,NPT
   38   V(I)= Y(I)
      END IF
C                                                *****************
C                                                **  LOG. SCALE **
C                                                *****************
      IF(IABS(IL).EQ.2) THEN
         DO 40 I=1,NPT
         IF(V(I).GT.0.) THEN
            V(I)= ALOG10(V(I))
         ELSE
            V(I)= -10.
         END IF
   40    CONTINUE
      END IF
C                                **************************************
C                                ** SET A NEW SCALE AND DRAW A FRAME.**
C                                **************************************
      IF(IPLOT.EQ.2) THEN
         ymax= -1.e10
         ymin=  1.e10
c
         do 50 i= 1,npt
         ymax= amax1(ymax,v(i))
         ymin= amin1(ymin,v(i))
   50    continue
c
         if(ymin.ge.0.) then
           ymax= 1.1*ymax
           ymin= 0.
         else
           ymax= amax1(0.,ymax)
           ymin= 1.1*ymin
         end if
      END IF
C
      IF(YMAX.LE.YMIN) YMAX= YMIN+1.0
      IF(IABS(IL).EQ.2) THEN
         IF(YMAX.GT.0.0) YMAX= YMAX+1.0
      END IF
C
      DX= (XMAX-XMIN)/XCM(I1)
      DY= (YMAX-YMIN)/YCM(J1)
      X0= XMIN
      Y0= YMIN
C
      CALL SCALEX (PL(I1),QL(J1),X0,Y0,DX,DY,ISC)
C
      PL1(ISC)= PL(I1)
      PR1(ISC)= PR(I1)
      QL1(ISC)= QL(J1)
      QR1(ISC)= QR(J1)
      XMIN1(ISC)= XMIN
      XMAX1(ISC)= XMAX
      YMAX1(ISC)= YMAX
      YMIN1(ISC)= YMIN
C                                                      *************
C                                                      **  FRAME. **
C                                                      *************
      CALL PLOT (PL(I1),QL(J1),3)
      CALL PLOT (PL(I1),QR(J1),2)
      CALL PLOT (PR(I1),QR(J1),2)
      CALL PLOT (PR(I1),QL(J1),2)
      CALL PLOT (PL(I1),QL(J1),2)
C                                                    ******************
C                                                    **  TICK MARKS. **
C                                                    ******************
      SCX= XCM(I1)/(nxtick+1)
      SCY= YCM(J1)/(nytick+1)
C
      X0= PL(I1)
      Y1= QL(J1)
      Y4= QR(J1)
      Y2= Y1 +0.25
      Y3= Y4 -0.25
C
      DO 62 K=1,nxtick
      X0= X0 +SCX
      CALL PLOT (X0,Y1,3)
      CALL PLOT (X0,Y2,2)
      CALL PLOT (X0,Y3,3)
      CALL PLOT (X0,Y4,2)
   62 CONTINUE
C
      Y0= QL(J1)
      X1= PL(I1)
      X4= PR(I1)
      X2= X1 +0.25
      X3= X4 -0.25
C
      DO 63 K=1,nytick
      Y0= Y0 +SCY
      CALL PLOT (X1,Y0,3)
      CALL PLOT (X2,Y0,2)
      CALL PLOT (X3,Y0,3)
      CALL PLOT (X4,Y0,2)
   63 CONTINUE
C                                                     **************
C                                                     ** NUMBERS. **
C                                                     **************
C
      hhs= 0.6
      CALL NUMBER (PL(I1)-0.8,QL(J1)-0.45,hhs,XMIN,0.,101)
      CALL NUMBER (PR(I1)-1.1,QL(J1)-0.45,hhs,XMAX,0.,101)
C
      CALL NUMBER (PL(I1)-1.8,QL(J1)     ,hhs,YMIN,0.,101)
      CALL NUMBER (PL(I1)-1.8,QR(J1)-0.30,hhs,YMAX,0.,101)
C
C                                                     **************
C                                                     **  LABELS. **
C                                                     **************
      XC= 0.5*(PL(I1)+PR(I1))
      XU= XC -1.60
      XD= XC -0.20*N2/2
C
      YR= QR(J1)+0.15
      YL= QL(J1)-0.70
C
      CALL SYMBOL (XU,YR,HH,LAB1,0.,N1)
      CALL SYMBOL (XD,YL,HH,LAB2,0.,N2)
C
      XL= PL(I1)-1.50
      YC= 0.5*(QL(J1)+QR(J1))
      CALL SYMBOL (XL,YC,HH,LAB3,0.,N3)
C                                     **********************************
C                                     **  NO PLOT IS MADE IF NPT1 < 0 **
C                                     **********************************
   70 IF(NPT1.LT.0) RETURN
C
      isk= 0
      if(iskip.ne.0) then
        call set_width (2.0)
      end if
c
      if(iskip.eq.1) then
        itot= 5
        ibr = 3 ! 2
      end if
      if(iskip.eq.2) then
        itot= 13
        ibr = 9
      end if
      if(iskip.eq.7) then
        itot= 3
        ibr = 1
      end if
c
      CALL PLOTL (U(1),V(1),ISC,3)
C**
      IF(IPLOT.EQ.1) THEN
         DO 100 I=1,NPT
         isk= isk +1
c
         if(iskip.eq.0) then
           CALL PLOTL (U(I),V(I),ISC,2)
         else
           if(mod(isk,itot).lt.ibr) then
             CALL PLOTL (U(I),V(I),ISC,2)
           else
             CALL PLOTL (U(I),V(I),ISC,3)
           end if
         end if
  100    CONTINUE
      ELSE
         DO 120 I=1,NPT-1
         CALL PLOTL (U(I+1),V(I)  ,ISC,2)
         CALL PLOTL (U(I+1),V(I+1),ISC,2)
  120    CONTINUE
      END IF
C**
      CALL PLOTL (U(NPT),V(NPT),ISC,3)
      call set_width (1.0)
C
      RETURN
      END
c
c
c------------------------------------
       subroutine set_width (awid)
c------------------------------------
       write(77,*) 'stroke'
       write(77,10) awid
   10  format(f4.1,' setlinewidth')
       return
       end
C
C
C-----------------------------------------------------------------------
      SUBROUTINE CPLOT3 (Q,xmax8,ymax8,zmax8,CHAR,NC)
C-----------------------------------------------------------------------
C***********************************************************************
C*   CONTOUR PLOTS OF SCALAR QUANTITIES.                               *
C***********************************************************************
      include    'param-spinRL.h'
      parameter  (mx1=meshx+1,my1=meshy+1,mz1=meshz+1)
      parameter  (mx=meshx,my=meshy,mz=meshz)
C
      CHARACTER*8  CHAR,LABEL
      COMMON/HEADR1/ LABEL,CDATE
      COMMON/HEADR2/ time
      integer*4      pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(mx1),pxc(mx1),pxl(mx1),pyr(my1),pyc(my1),
     *               pyl(my1),pzr(mz1),pzc(mz1),pzl(mz1)
      real*8         Q(mx,my,mz),xmax8,ymax8,zmax8
      real*4         a(2048),b(2048),WW(2048),cut(200,4)
C
      J0= my/2 +1
      K0= mz/2 +1
      xmax= xmax8
      ymax= ymax8
      zmax= zmax8
C
C* 1. PLOT AT K= K0: SUBSCRIPT J FIRST.
C
      npx= 0
      IJ= 0
      qc = 1./16.
c***
      do 10 i= 1,mx
      npx= npx +1
      IR= PXR(I)
      IL= PXL(I)
c
      npy= 0
      DO 10 j= 1,my
      npy= npy +1
      JR= PYR(J)
      JL= PYL(J)
C
      IJ= IJ+1
      a(IJ)= qc*(   Q(IR,JR,K0) +2.*Q(IR,J,K0)    +Q(IR,JL,K0)
     *          +2.*Q(I ,JR,K0) +4.*Q(I ,J,K0) +2.*Q(I ,JL,K0)
     *          +   Q(IL,JR,K0) +2.*Q(IL,J,K0)    +Q(IL,JL,K0) )
   10 CONTINUE
C
C* 2. PLOT AT J= J0: SUBSCRIPT K FIRST.
C
      npx= 0
      IJ= 0
      qc = 1./16.
c***
      do 20 i= 1,mx
      npx= npx +1
      IR= PXR(I)
      IL= PXL(I)
c
      npz= 0
      DO 20 k= 1,mz
      npz= npz +1
      KR= PZR(K)
      KL= PZL(K)
C
      IJ= IJ+1
      b(IJ)= qc*(   Q(IR,J0,KR) +2.*Q(IR,J0,K)    +Q(IR,J0,KL)
     *          +2.*Q(I ,J0,KR) +4.*Q(I ,J0,K) +2.*Q(I ,J0,KL)
     *          +   Q(IL,J0,KR) +2.*Q(IL,J0,K)    +Q(IL,J0,KL) )
   20 CONTINUE
C
C
      HH = 0.70
      CALL SYMBOL (0.1,18.2,HH,LABEL,0.,8)
      CALL SYMBOL (13.0,0.7,HH,CDATE,0.,9)
      CALL SYMBOL (13.0,0.1,HH,'T =',0.,3)
      CALL VALUES (999.0,999.0,HH,TIME,0.,101)
C
      XL1=  1.8
      XR1=  9.3
      XL2= 10.0
      XR2= 17.5
C
      ZL=  1.0
      ZR=  ZL +(XR1 -XL1)*ZMAX/XMAX
      if(zr.gt.10.) zr= 10.
c                          <--- Limit elongated y-length.
C
      YL=  ZR +1.
      YR=  YL +(XR1 -XL1)*YMAX/XMAX
      if(yr.gt.25.) yr= 25.
c                          <--- Limit elongated y-length.
C
      XC1= 0.5*(XR1+XL1)
      XC2= 0.5*(XR2+XL2)
      YC=  0.5*(YR+YL)
      ZC=  0.5*(ZR+ZL)
C
C---------------------------------------------
C*  **MAXIMUM OF THE VECTORS**
C---------------------------------------------
C
      AM2= 0.
      AM4= 0.
C
      DO 100 IJ= 1,npx*npy
      AM2= AMAX1(AM2,ABS(A(IJ)))
  100 CONTINUE
C
      DO 200 IJ= 1,npx*npz
      AM4= AMAX1(AM4,ABS(B(IJ)))
  200 CONTINUE
C
      AMS= AMAX1(AM2,AM4)
      IF(AMS.LT.1.E-10) AMS=999.0
C
      CALL SYMBOL (ZL,0.10,HH,'SCALAR.MAX= ',0.,12)
      CALL VALUES (999.0,999.0,HH,AMS,0.,101)
C
C---------------------------------------------
C*  (1): Contours in (x,z) plane.
C---------------------------------------------
C
      CALL SETSCL (0.,0.,ZMAX,XMAX,ZL,XL2,ZR,XR2,GDZ,GDX,
     &             NC,CHAR,6,' (X-Z)',0.4,
     &             1,'Z',0.4,1,'X',0.4,1)
C
      CALL VALUES (ZL-0.45,XL2-0.5,HH,0.0,0.,101)
      CALL VALUES (ZL-1.3,XR2-0.3, HH,XMAX,0.,101)
      CALL VALUES (ZR-1.3,XL2-0.5, HH,ZMAX,0.,101)
C
      NXZ= npx*npz
      CALL DAISHO (B,NXZ,WAMIN,WAMAX)
C
      NCONTR= 11
      CALL EQCNTR (B,WW,npz,npx,ZL,XL2,ZR,XR2,WAMIN,0.0,WAMAX,
     *             NCONTR,1)
C
C---------------------------------------------
C*  (2): Contours in (x,y) plane.
C---------------------------------------------
C
      CALL SETSCL (0.,0.,YMAX,XMAX,YL,XL2,YR,XR2,GDY,GDX,
     &             1,' ',1,' ',0.4,
     &             1,'Y',0.4,1,'X',0.4,1)
C
      CALL VALUES (YL-0.45,XL2-0.5,HH,0.0,0.,101)
      CALL VALUES (YL-1.3,XR2-0.3, HH,XMAX,0.,101)
      CALL VALUES (YR-0.3,XL2-0.5, HH,YMAX,0.,101)
C
      NXY= npx*npy
      CALL DAISHO (A,NXY,WAMIN,WAMAX)
C
      NCONTR= 11
      CALL EQCNTR (A,WW,npy,npx,YL,XL2,YR,XR2,WAMIN,0.0,WAMAX,
     *             NCONTR,1)
C
C---------------------------------------------
C*  (3): Cut plots.
C---------------------------------------------
C
      do 300 jj= 1,4
      j= (my/4)*(jj-1) +1
c
      do 300 i= 1,mx
      cut(i,jj)= q(i,j,k0)
  300 continue
c
c
      amax7= -1.e+10
      amin7=  1.e+10
c
      do 320 jj= 1,4
      do 320 i= 1,mx
      amax7= amax1(cut(i,jj),amax7)
      amin7= amin1(cut(i,jj),amin7)
  320 continue
c
      if(amax7.lt.0.) amax7= 0.
      if(amin7.gt.0.) amin7= 0.
c
c
      dd= amax7 -amin7
      dx= (YR -YL)/6.
      dy= XR1 -XL1
c
      do 340 jj= 1,4
      xo= YL +1.5*(jj-1)*dx
      xu= xo +dx
      call plot (xo,XL1,3)
      call plot (xu,XL1,2)
      call plot (xu,XR1,2)
      call plot (xo,XR1,2)
      call plot (xo,XL1,2)
c
c* zero line.
      x1= xo +dx*(0. -amin7)/dd
      call plot (x1,XL1,3)
      call plot (x1,XR1,2)
c
      x1= xo +dx*(cut(1,jj) -amin7)/dd
      y1= XL1 
      call plot (x1,y1,3)
c
      do 340 i= 1,mx
      x1= xo  +dx*(cut(i,jj) -amin7)/dd
      y1= XL1 +dy*(i-1)/float(mx-1)
c
      call plot (x1,y1,2)
  340 continue
c---------------------
      call chart
c---------------------
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE EQCNTR(U,W,NX,NY,XL,YL,XR,YR,UMIN,UBUND,UMAX,
     &                  LANK,IWAKU)
C-----------------------------------------------------------------------
C  << EQCNTR >>
C         PRESENTED BY KUNIHIKO.WATANABE 14.NOV.1989
C         REVICED   BY HISANORI.TAKAMARU 16.MAR.1990
C-----------------------------------------------------------------------
C     1. FUNCTION
C        (1) TO DRAW TOKOSEN
C     2. ARGUMENTS   (SIZE)   (I/O)     (MEANING)
C        (1) U       NX,NY     (I)       WORLD VALUE
C        (2) W       NX,NY     (I)       WORK ARRAY (REAL*4)
C        (3) XL,XR,YL,YR       (I)       ABSOLUTE COORDINATE VALUE
C        (4) UMIN,UMAX         (I)       HEIGHT OF MAX & MIN
C                                        UMIN>UMAX : AUTOMATIC CONTROL
C        (5) UBUND             (I)       DRAW DASH LINE (U < UBUND)
C        (6) LANK              (I)       NUMBER OF DRAW LINES
C        (7) IWAKU             (I)       =1 : DRAW FRAME
C     3. CALLED BY
C             (** NOTHING **)
C     4. CALLS
C             (** PLOT   **)
C-----------------------------------------------------------------------
      DIMENSION U(1),W(1)
C
      IF (NX.LT.2) RETURN
      IF (NY.LT.2) RETURN
      IF (XR.LT.XL) RETURN
      IF (YR.LT.YL) RETURN
C
      NXY = NX*NY
      NXM1 = NX - 1
      NYM1 = NY - 1
C
      DX = (XR-XL)/ FLOAT(NXM1)
      DY = (YR-YL)/ FLOAT(NYM1)
C
      UMAX1 = UMAX
      UMIN1 = UMIN
C
      IF(UMAX1.GT.(1.000001*UMIN1)) THEN
C
        DO 10 I = 1 , NXY
          W(I) = U(I) - UMIN1
          IF(U(I).GT.UMAX1) W(I) = UMAX1 - UMIN1
          IF(U(I).LT.UMIN1) W(I) = 0.
   10   CONTINUE
C
      ELSE
C
        UMAX1=-1.E+30
        UMIN1= 1.E+30
        DO 20 I = 1 , NXY
          UMAX1=AMAX1(UMAX1,U(I))
          UMIN1=AMIN1(UMIN1,U(I))
   20   CONTINUE
        DO 25 I = 1 , NXY
          W(I) = U(I) - UMIN1
   25   CONTINUE
C
      ENDIF
C
C------------------------------------------------
      IF(UMAX1.LE.(1.000001*UMIN1))  RETURN
C------------------------------------------------
C
      IF(IWAKU.EQ.1) THEN
        CALL PLOT(XL,YL,3)
        CALL PLOT(XR,YL,2)
        CALL PLOT(XR,YR,2)
        CALL PLOT(XL,YR,2)
        CALL PLOT(XL,YL,2)
        CALL PLOT(XL,YL,3)
      ENDIF
C
      ULD = FLOAT(LANK+1) / (UMAX1-UMIN1)
      EPS = 1.0E-8
C
      NXYM1 = NXM1*NYM1
      DO 9000  IJNXY1 = 1,NXYM1
        J = (IJNXY1-1)/NXM1 + 1
        I = IJNXY1 - (J-1)*NXM1
C
          I1 = I + NX * (J - 1)
          I2 = I1 + 1
          I3 = I1 + 1 + NX
          I4 = I1 + NX
C
          U1 =  W(I1) * ULD
          U2 =  W(I2) * ULD
          U3 =  W(I3) * ULD
          U4 =  W(I4) * ULD
C
          K1 = IFIX(U1)
          K2 = IFIX(U2)
          K3 = IFIX(U3)
          K4 = IFIX(U4)
C
          J1 = IABS(K2-K1)
          J2 = IABS(K3-K2)
          J3 = IABS(K4-K3)
C
          IF(J1.NE.0) THEN
            DO 1000 LL = 1 , J1
              U0 = FLOAT(LL) + FLOAT(MIN0(K1,K2))
                UJOUGE = U0/ULD + UMIN1
                IF (UJOUGE.LT.UBUND) THEN
                  JOUGE = 4
                ELSE
                  JOUGE = 1
                END IF
C
              IF(ABS(U2-U1).LT.EPS)                 GO TO 1000
C
              X1 = XL + DX * ( (U0-U1)/(U2-U1) + FLOAT(I-1) )
              Y1 = YL + DY * FLOAT(J-1)
C
              IF( ((U3-U0)*(U2-U0)).GT.0. )         GO TO 1100
              IF( ( (U0-U2).GT.0. ).AND.( (U0-U4).GT.0. ) ) GO TO 1100
              IF( ABS(U3-U2).LT.EPS )               GO TO 1100
C
                X2 = XL + DX * FLOAT(I)
                Y2 = YL + DY * ( (U0-U2)/(U3-U2) + FLOAT(J-1) )
C
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
C
 1100         CONTINUE
              IF( ((U4-U0)*(U3-U0)).GT.0. )         GO TO 1200
              IF( ((U1-U0)*(U3-U0)).GT.0. )         GO TO 1200
              IF( ((U2-U0)*(U4-U0)).GT.0. )         GO TO 1200
              IF( ABS(U4-U3).LT.EPS )               GO TO 1200
C
                X2 = XL + DX * ( (U0-U4)/(U3-U4) + FLOAT(I-1) )
                Y2 = YL + DY * FLOAT(J)
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
C
 1200         CONTINUE
              IF( ((U1-U0)*(U4-U0)).GT.0. )         GO TO 1300
              IF( ( (U0-U1).GT.0. ).AND.( (U0-U3).GT.0. ) ) GO TO 1300
              IF( ABS(U1-U4).LT.EPS )               GO TO 1300
C
                X2 = XL + DX * FLOAT(I-1)
                Y2 = YL + DY*((U0-U1)/(U4-U1)+FLOAT(J-1))
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
 1300         CONTINUE
 1000       CONTINUE
C
          ENDIF
C
          IF(J2.NE.0) THEN
C
            DO 2000 LL = 1 , J2
              U0 = FLOAT(LL) + FLOAT(MIN0(K2,K3))
                UJOUGE = U0/ULD + UMIN1
                IF (UJOUGE.LT.UBUND) THEN
                  JOUGE = 4
                ELSE
                  JOUGE = 1
                END IF
              IF( ABS(U3-U2).LT.EPS )               GO TO 2000
C
              X1 = XL + DX * FLOAT(I)
              Y1 = YL + DY * ( (U0-U2)/(U3-U2) + FLOAT(J-1) )
C
              IF( ((U4-U0)*(U3-U0)).GT.0. )         GO TO 2100
              IF( ( (U0-U1).GT.0. ).AND.( (U0-U3).GT.0. ) ) GO TO 2100
              IF( ABS(U4-U3).LT.EPS )               GO TO 2100
C
                X2 = XL + DX * ( (U0-U4)/(U3-U4) + FLOAT(I-1) )
                Y2 = YL + DY * FLOAT(J)
C
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
C
 2100         CONTINUE
              IF( ((U1-U0)*(U4-U0)).GT.0. )         GO TO 2200
              IF( ((U1-U0)*(U3-U0)).GT.0. )         GO TO 2200
              IF( ((U2-U0)*(U4-U0)).GT.0. )         GO TO 2200
              IF( ABS(U1-U4).LT.EPS )               GO TO 2200
C
                X2 = XL + DX * FLOAT(I-1)
                Y2 = YL + DY * ( (U0-U1)/(U4-U1)+FLOAT(J-1) )
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
 2200         CONTINUE
 2000       CONTINUE
C
          ENDIF
C
          IF(J3.NE.0) THEN
C
            DO 3000 LL = 1 , J3
              U0 = FLOAT(LL) + FLOAT(MIN0(K3,K4))
                UJOUGE = U0/ULD + UMIN1
                IF (UJOUGE.LT.UBUND) THEN
                  JOUGE = 4
                ELSE
                  JOUGE = 1
                END IF
              IF( ABS(U4-U3).LT.EPS )               GO TO 3000
C
              X1 = XL + DX * ( (U0-U4)/(U3-U4) + FLOAT(I-1) )
              Y1 = YL + DY * FLOAT(J)
C
              IF( ((U1-U0)*(U4-U0)).GT.0. )         GO TO 3100
              IF( ( (U0-U2).GT.0. ).AND.( (U0-U4).GT.0. ) ) GO TO 3100
              IF( ABS(U1-U4).LT.EPS )               GO TO 3100
C
                X2 = XL + DX * FLOAT(I-1)
                Y2 = YL + DY * ( (U0-U1)/(U4-U1) + FLOAT(J-1) )
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
 3100         CONTINUE
 3000       CONTINUE
          ENDIF
 9000 CONTINUE
C
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE SETSCL (WMINX,WMINY,WMAXX,WMAXY, XL,YL,XR,YR, GDX,GDY,
     &           N1,CHAR1, N2,CHAR2, HIGHT1,
     &           NNX,CHARX, HIGHTX, NNY,CHARY, HIGHTY, IWAKU)
C-----------------------------------------------------------------------
C  << SETSCL >>                   /CHAR1/
C                          WMAXY  +--------------------+  (XL,YL)
C                               Y |             (XR,YR)|  (XR,YR) ON 0
C                               R |                    |
C                               A |                    |
C    (WMINX,WMINY)              H |                    |
C    (WMAXX,WMAXY) ON IS        C |                    |
C                                 |(XL,YL)             |
C                          WMINY  +--------+--+--------+
C                                 WMINX  /CHARX/       WMAXX
C-----------------------------------------------------------------------
*
C     SETSCL
C
C     1. FUNCTION
C        (1) TO SCALE THE GRAPHICS BY CALCOMP SPECIFICATIONS
C     2. ARGUMENTS            (I/O)     (MEANING)
C        (1) WMINX,WMAXX,
C            WMINY,WMAXY       (I)       WORLD COORDINATE VALUE
C        (2) XL,XR,YL,YR       (I)       ABSOLUTE COORDINATE VALUE
C        (3) GDX,GDY           (O)       SCALING FACTOR OF COORDINATE
C                                        FROM WORLD TO ABSOLUTE
C        (4) CHAR1,CHARX,CAHRY (I)       TITLE ON GRAPH,X-AXIS,Y-AXIS
C        (5) IWAKU             (I)       DRAW FRAME (0:OFF ; 1:ON)
C                                         999 : WRITE OUT ONLY TITLE,
C                                                    NOT DRAW OTHERWISE
C     3. CALLED BY
C             (** NOTHING **)
C     4. CALLS
C             (** PLOT   **)
C             (** SYMBOL **)
C             (** NUMBER **)
C-----------------------------------------------------------------------
      CHARACTER*1  CHAR1(1),CHAR2(1),CHARX(1),CHARY(1)
*
      IF (WMAXX.LE.WMINY) GOTO 9999
      IF (WMAXX.LE.WMINY) GOTO 9999
      IF (XR.LE.XL)       GOTO 9999
      IF (YR.LE.YL)       GOTO 9999
*
      GDX= (XR-XL)/(WMAXX-WMINX)
      GDY= (YR-YL)/(WMAXY-WMINY)
*
      XC = 0.5*( XR + XL )
      YC = 0.5*( YR + YL )
*
      IF (N1 .GT.0) THEN
        IF (HIGHT1.GT.0) THEN
          XS1= XC -0.5*N1*HIGHT1
          XS2= XS1 +(N1+1)*HIGHT1
          CALL SYMBOL(XS1,YR+0.1,HIGHT1,CHAR1(1),0.,N1)
          CALL SYMBOL(XS2,YR+0.1,HIGHT1,CHAR2(1),0.,N2)
        END IF
      END IF
C-----------------------------------------------------------------------
      IF (IWAKU.EQ.999) RETURN
C-----------------------------------------------------------------------
*
      IF (IWAKU.EQ.1) THEN
        CALL PLOT (XL,YL,3)
        CALL PLOT (XL,YR,2)
        CALL PLOT (XR,YR,2)
        CALL PLOT (XR,YL,2)
        CALL PLOT (XL,YL,2)
        CALL PLOT (999.,999.0,3)
      END IF
*
      IF (NNX.GT.0) THEN
        IF (HIGHTX.GT.0) THEN
          CALL SYMBOL(XC-0.5*HIGHTX*NNX,YL-0.5,HIGHTX,CHARX(1),0.,1)
          DO 200 NNX1=2,NNX
  200     CALL SYMBOL(999.0,999.0,HIGHTX,CHARX(NNX1),0.,1)
        END IF
      END IF
      IF (NNY.GT.0) THEN
        IF (HIGHTY.GT.0) THEN
          CALL SYMBOL(XL-0.5,YC-0.5*HIGHTY*NNY,HIGHTY,CHARY(1),0.,1)
          DO 300 NNY1=2,NNY
  300     CALL SYMBOL(999.0,999.0,HIGHTY,CHARY(NNY1),0.,1)
        END IF
      ELSE IF(NNY.LT.0) THEN
        IF (HIGHTY.GT.0) THEN
          CALL SYMBOL(XC-0.5*HIGHTY*NNY,YC,HIGHTY,CHARY(1),0.,1)
          DO 400 NNY1=2,NNY
  400     CALL SYMBOL(999.0,999.0,HIGHTY,CHARY(NNY1),0.,1)
        END IF
      END IF
*
      RETURN
*
C-----------------------------------------------------------------------
*
 9999 CONTINUE
      WRITE(6,*) '**********  ABNORMAL WORLD COORDINATE ********'
      WRITE(6,*) '      '
      WRITE(6,*) '    WMAXX =',WMAXX,' WMINX = ',WMINX
      WRITE(6,*) '    WMAXY =',WMAXY,' WMINY = ',WMINY
      WRITE(6,*) '    XL,YL,XR,YR =',XL,YL,XR,YR
      WRITE(6,*) '    FCTR  =',FCTR
      WRITE(6,*) '      '
      CALL CHART
      CALL SYMBOL(1.0,10.0,0.2,' ABNORMAL WORLD COORDINATE CALL',0.,31)
      CALL SYMBOL(1.0,09.0,0.2,' WMAXX =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,WMAXX,0.,2)
      CALL SYMBOL(1.0,08.5,0.2,' WMINX =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,WMINY,0.,2)
      CALL SYMBOL(1.0,08.0,0.2,' WMAXY =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,WMAXY,0.,2)
      CALL SYMBOL(1.0,07.5,0.2,' WMINY =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,WMINY,0.,2)
      CALL SYMBOL(1.0,07.0,0.2,' FCTR  =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,FCTR,0.,2)
      CALL SYMBOL(1.0,06.5,0.2,' XLEFT =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,XL,0.,2)
      CALL SYMBOL(1.0,06.0,0.2,' YLEFT =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,YL,0.,2)
      CALL SYMBOL(1.0,05.5,0.2,' XRIGHT=',0.,8)
      CALL NUMBER(999.0,999.0,0.2,XR,0.,2)
      CALL SYMBOL(1.0,05.0,0.2,' YRIGHT=',0.,8)
      CALL NUMBER(999.0,999.0,0.2,YR,0.,2)
      CALL EXIT
      RETURN
      END
c
c
C-----------------------------------------------------------------------
      SUBROUTINE SCALEX (XCM,YCM,X00,Y00,DX,DY,ISC)
C-----------------------------------------------------------------------
      COMMON/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
C
      X0(ISC)= X00
      Y0(ISC)= Y00
      DXI(ISC)= 1./DX
      DYI(ISC)= 1./DY
C
      XL(ISC)= XCM
      YL(ISC)= YCM
C
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE PLOTL (X,Y,ISC,IPL)
C-----------------------------------------------------------------------
      COMMON/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
C
      XCM= XL(ISC) +DXI(ISC)*(X -X0(ISC))
      YCM= YL(ISC) +DYI(ISC)*(Y -Y0(ISC))
C
      CALL PLOT (XCM,YCM,IPL)
C
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE VALUES (X,Y,HEIGHT,VAL,THETA,IFMAT)
C-----------------------------------------------------------------------
C  << VALUES >>
C     1. FUNCTION
C        (1) TO DRAW VARIABLE
C     2. ARGUMENTS   (SIZE)   (I/O)     (MEANING)
C        (1) X,Y               (I)       ABSOLUTE COORDINATE VALUE
C        (2) HEIGHT            (I)       DRAW OUT SIZE ON PAPER
C        (3) VAL               (I)       VARIABLE
C        (4) THETA             (I)       ANGLE
C        (5) IFMAT             (I)       FORMAT TYPE
C     3. CALLED BY
C             (** NOTHING **)
C     4. CALLS
C             (** NUMBER **)
C             (** SYMBOL **)
C-----------------------------------------------------------------------
*        IFMAT = (N100)*100 + KETA
*        N100 = 0 : INTEGER FORMAT
*        N100 = 1 : F FORMAT ::  NUMBER(X,Y,HEIGHT,VAL,THETA,KETA)
*        N100 = 2 : E FORMAT ::
*        N100 = 3 : POWER OF TEN FORMAT
*        N100 = OTHEWISE : NOT WRITE OUT
*-----------------------------------------------------------------------
*
      REAL*4 VAL
      CHARACTER CHR13*13,CHR12*12,CHR3*3
      CHARACTER*1 MINUS,ZERO,BLANK
      PARAMETER(RATIO = 6./7. )
      DATA MINUS/'-'/,ZERO/'0'/,BLANK/' '/
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (IFMAT.LT.0) RETURN
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      N100 = IFMAT/100
      KETA = IFMAT - N100*100
*
      IF (N100.EQ.0) THEN
        CALL NUMBER(X,Y,HEIGHT,VAL,THETA,-1)
      ELSE IF (N100.EQ.1) THEN
        CALL NUMBER(X,Y,HEIGHT,VAL,THETA,KETA)
      ELSE IF (N100.EQ.2) THEN
        CHR13 = '             '
        CHR12 = '            '
        IF (KETA.EQ.0) THEN
          WRITE(CHR13,'(1PE13.6)') VAL
          CHR12(1:4) = CHR13(1:3)//'E'
          NUMSYM = 4
        ELSE
          KETA = KETA + 1
          IF (VAL.LT.0.) THEN
            CHRVAL = VAL - 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+3) = CHR13(1:KETA+2)//'E'
            NUMSYM = KETA + 3
          ELSE IF (VAL.EQ.0) THEN
            CHRVAL = VAL
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+3) = CHR13(1:KETA+2)//'E'
            NUMSYM = KETA + 3
          ELSE
            CHRVAL = VAL + 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+2) = CHR13(2:KETA+2)//'E'
            NUMSYM = KETA + 2
          END IF
        END IF
        CHR3 = '   '
*
        IF (CHR13(11:11) .EQ. MINUS) THEN
          IF (CHR13(12:12) .EQ. ZERO  .OR.
     &        CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:2) = '-'//CHR13(13:13)
          ELSE
            CHR3(1:3) = '-'//CHR13(12:13)
          END IF
          NUMSY1 = 3
        ELSE
          IF (CHR13(12:12) .EQ. ZERO  .OR.
     &        CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:1) = CHR13(13:13)
            NUMSY1 = 1
          ELSE
            CHR3(1:2) = CHR13(12:13)
            NUMSY1 = 2
          END IF
        END IF
        AKAKU = 2. * 3.1415927 / 360.
        COST = COS(THETA*AKAKU)
        CALL SYMBOL(X,Y,HEIGHT,CHR12,THETA,NUMSYM)
        CALL SYMBOL(999.,999.,HEIGHT,CHR3,THETA,NUMSY1)
      ELSE IF (N100.EQ.3) THEN
        CHR13 = '             '
        CHR12 = '            '
        IF (KETA.EQ.0) THEN
          WRITE(CHR13,'(1PE13.6)') VAL
          CHR12(1:6) = CHR13(1:3)//'X10'
          NUMSYM = 6
        ELSE
          KETA = KETA + 1
          IF (VAL.LT.0.) THEN
            CHRVAL = VAL - 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+5) = CHR13(1:KETA+2)//'X10'
            NUMSYM = KETA + 5
          ELSE
            CHRVAL = VAL + 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+4) = CHR13(2:KETA+2)//'X10'
            NUMSYM = KETA + 4
          END IF
        END IF
        CHR3 = '   '
*
        IF (CHR13(11:11) .EQ. MINUS) THEN
          IF (CHR13(12:12) .EQ. ZERO  .OR.
     &        CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:2) = '-'//CHR13(13:13)
          ELSE
            CHR3(1:3) = '-'//CHR13(12:13)
          END IF
          NUMSY1 = 3
        ELSE
          IF (CHR13(12:12) .EQ. ZERO  .OR.
     &        CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:1) = CHR13(13:13)
            NUMSY1 = 1
          ELSE
            CHR3(1:2) = CHR13(12:13)
            NUMSY1 = 2
          END IF
        END IF
        AKAKU = 2. * 3.1415927 / 360.
        COST = COS(THETA*AKAKU)
        SINT = SIN(THETA*AKAKU)
        CALL SYMBOL(X,Y,HEIGHT,CHR12,THETA,NUMSYM)
*
*                                             *******************
*                                             ** EXPONENT PART **
*                                             *******************
*
        H2 = HEIGHT * 5./7.
        X1 = (NUMSYM+1)* HEIGHT * RATIO
        Y1 = HEIGHT * 4./7.
        IF (ABS(THETA).LT.1E-04) THEN
          X1 = X + X1
          Y1 = Y + Y1
        ELSE
          X2 =     X1 * COST - Y1 * SINT
          Y1 = Y + X1 * SINT + Y1 * COST + H2*COST
          X1 = X + X2                    - H2*SINT
        END IF
        CALL SYMBOL(X1,Y1,H2,CHR3,THETA,NUMSY1)
      END IF
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE WDASH (X1,Y1,X2,Y2,IPEN )
C-----------------------------------------------------------------------
*  << WDASH  >>                      VER 2.00   16.MAR.1990
C
C     1. FUNCTION
C        (1) TO DRAW LINE FROM (X1,Y1) TO (X2,Y2) BY WDASH
C                            IN ABSOLUTE COORDINATE
C     2. ARGUMENTS            (I/O)     (MEANING)
C        (1) X1,X2,Y1,Y2       (I)       ABSOLUTE COORDINATE VALUE
C        (2) IPEN              (I)       PEN TYPE OF 'WDASH'
C     3. CALLED BY
C             (** EQCNTR  **)
C             (** WDASHL  **)
C     4. CALLS
C             (** PLOT   **)
C-----------------------------------------------------------------------
C       IPEN : MEANING           - : 0.05 (CM)
C        1   :       LINE     -------------------
C        2   :  DASH LINE     --- --- --- --- ---
C        3   :  DASH LINE     -- -- -- -- -- -- --
C        4   :  DASH LINE     - - - - - - - - - -
C        5   :  1 POINT DASH  ---- - ---- - ---- -
C        6   :  2 POINT DASH  --2.0-- - - --2.0--
C   OTHERWISE:  LINE          ---------------------
C-----------------------------------------------------------------------
C
      H1  =  0.05
      H2  =  2.0 * H1
      H3  =  3.0 * H1
      H4  =  4.0 * H1
      H20 = 20.0 * H1
      CALL PLOT ( X1 , Y1 , 3 )
      K = - 1
      IF(IPEN.LT.2) THEN
        GO TO 999
      ELSE IF(IPEN.EQ.2) THEN
        HH1 = H3
        HH2 = H1
      ELSE IF (IPEN.EQ.3) THEN
        HH1 = H2
        HH2 = H1
      ELSE IF (IPEN.EQ.4) THEN
        HH1 = H1
        HH2 = H1
      ELSE IF (IPEN.EQ.5) THEN
        HH1 = H4
        HH2 = H1
        HH3 = H1
        HH4 = H1
      ELSE IF (IPEN.EQ.6) THEN
        HH1 = H20
        HH2 = H1
        HH3 = H1
        HH4 = H1
        HH5 = H1
        HH6 = H1
      END IF
      IF(IPEN.LT.5) THEN
        RLENG = SQRT ( ( X2 - X1 ) **2 + ( Y2 - Y1 ) **2 )
        IF(RLENG.LT.1.0E-5) GOTO 999
        IF(RLENG.LT.HH1) GOTO 999
        COSTH = ( X2 - X1 ) / RLENG
        SINTH = ( Y2 - Y1 ) / RLENG
        D = HH1
        X = X1 + D * COSTH
        Y = Y1 + D * SINTH
        CALL PLOT ( X , Y , ( 5 + K ) / 2 )
        K = - K
        D = D + HH2
        HHH = HH1
        HH1 = HH2
        HH2 = HHH
  200   IF(D.LE.RLENG) THEN
          X = X1 + D * COSTH
          Y = Y1 + D * SINTH
          CALL PLOT ( X , Y , ( 5 + K ) / 2 )
          K = - K
          HHH = HH1
          HH1 = HH2
          HH2 = HHH
          D=D+HH1
          GOTO 200
        END IF
      ELSE IF (IPEN.EQ.5) THEN
        RLENG = SQRT ( ( X2 - X1 ) **2 + ( Y2 - Y1 ) **2 )
        IF(RLENG.LT.1.0E-5) GOTO 999
        IF(RLENG.LT.HH1) GOTO 999
        COSTH = ( X2 - X1 ) / RLENG
        SINTH = ( Y2 - Y1 ) / RLENG
        D = HH1
        X = X1 + D * COSTH
        Y = Y1 + D * SINTH
        CALL PLOT ( X , Y , ( 5 + K ) / 2 )
        K = - K
        D = D + HH2
        HHH = HH1
        HH1 = HH2
        HH2 = HH3
        HH3 = HH4
        HH4 = HHH
  500   IF(D.LE.RLENG) THEN
          X = X1 + D * COSTH
          Y = Y1 + D * SINTH
          CALL PLOT ( X , Y , ( 5 + K ) / 2 )
          K = - K
          HHH = HH1
          HH1 = HH2
          HH2 = HH3
          HH3 = HH4
          HH4 = HHH
          D=D+HH1
          GOTO 500
        END IF
      ELSE IF (IPEN.EQ.6) THEN
        RLENG = SQRT ( ( X2 - X1 ) **2 + ( Y2 - Y1 ) **2 )
        IF(RLENG.LT.1.0E-5) GOTO 999
        IF(RLENG.LT.HH1) GOTO 999
        COSTH = ( X2 - X1 ) / RLENG
        SINTH = ( Y2 - Y1 ) / RLENG
        D = HH1
        X = X1 + D * COSTH
        Y = Y1 + D * SINTH
        CALL PLOT ( X , Y , ( 5 + K ) / 2 )
        K = - K
        D = D + HH2
        HHH = HH1
        HH1 = HH2
        HH2 = HH3
        HH3 = HH4
        HH4 = HH5
        HH5 = HH6
        HH6 = HHH
  600   IF(D.LE.RLENG) THEN
          X = X1 + D * COSTH
          Y = Y1 + D * SINTH
          CALL PLOT ( X , Y , ( 5 + K ) / 2 )
          K = - K
          HHH = HH1
          HH1 = HH2
          HH2 = HH3
          HH3 = HH4
          HH4 = HH5
          HH5 = HH6
          HH6 = HHH
          D=D+HH1
          GOTO 600
        END IF
      END IF
  999 CALL PLOT ( X2 , Y2 , ( 5 + K ) / 2 )
      CALL PLOT ( X2 , Y2 , 3)
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE DAISHO(X  ,NX,XMIN1,XMAX1)
C-----------------------------------------------------------------------
      DIMENSION X(1)
C
      XMAX1= X(1)
      XMIN1= X(1)
      DO 100 I=2,NX
      XMAX1= AMAX1(XMAX1,X(I) )
      XMIN1= AMIN1(XMIN1,X(I) )
  100 CONTINUE
      RETURN
      END
c
c
c***************************************************************
c*     This program package generates a UNIX postscript        *
c*     graphic file when called by calcomp-compatible          *
c*     /plot23.f/.                                             *
c***************************************************************
c----------------------------------------------------------
c      PostScript header by fortran
c        T. Ogino (Nagoya University) February 27, 1992
c      Modified to conform GSIPP commands
c        Motohiko Tanaka (NIFS)       November 23, 1993
c
c----------------------------------------------- 5/27/96 -------
c     This PS-Adobe-2.0 header allows us full paging features in
c     the Ghostview.  To scroll up the page (backward), click the 
c     page number and press two buttons of mouse simultaneously.
c
c     Consult: A.Saitou (Kyoto U.)  The definition of /@eop  
c    needs stroke for line drawings (not in the TeX header).
c---------------------------------------------------------------
       subroutine gopen (nframe)
c----------------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
c
c*  This is an Adobe-2.0 postscript file.
c
       write(77,10)
   10  format('%!PS-Adobe-2.0',/
     *        '%%Pages: (atend)',/
     *        '%%PageOrder: Ascend',/
     *        '%%EndComments',/
     *        '%%BeginDocument')
c
c%%%%%%%%%%%%%%%%%%% Procedure Defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     write(77,11) 
c  11 format('%%BoundingBox: 150. 400. 550. 600.')
c
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/
     *       '/m {moveto} bind def  % x y m -- move to position')
c
      write(77,23) 
   23 format('/tr {/Times-Roman findfont} bind def',/
     *       '/sf {scalefont} bind def',/
     *       '/se {setfont} bind def',/
     *       '/ro {rotate}  bind def',/
     *       '/tl {translate} bind def',/
     *       '/sc {scale} bind def')
c
      write(77,24) 
   24 format('/@bop          % @bop -- begin the a new page',/
     *       '{erasepage newpath initgraphics',/
     *       '/SaveImage save def',/
     *       '} bind def')
c
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/
     *       '{stroke showpage',/
     *       ' SaveImage restore',/
     *       '} bind def')
c
      write(77,26) 
   26 format('/@end          % @end -- done the whole shebang',/
     *       ' /end load def')
c
      write(77,27) 
   27 format('/dir 0 def')
c
      write(77,29) 
   29 format('/s             % string s -- show the string',/
     *       '{dir 1 eq',/
     *       ' {gsave currentpoint translate 90 rotate 0 0 moveto',/
     *       ' show grestore}',/
     *       ' {show} ifelse',/
     *       '} bind def')
c
      write(77,31)
   31 format('%%EndDocument',/
     *       '%%EndProlog',/
     *       '%%BeginSetup',/
     *       '/Resolution 300 def',/
     *       '/#copies 1 def',/
     *       '%%EndSetup')
c
c%%%%%%%%%%%%%%%%%%% End of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
c
c*  initiate the page one.
c
       nfrm = nframe
c
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%Page:',1x,i2,1x,i2)
c
       write(77,30) 
   30  format('%%BeginPageSetup',/
     *        '%%EndPageSetup',/
     '        '@bop')
c
c
c*  Set magnifying factor (GSIPP to Sun coordinate).
c   Rotate and translate to output on A4-L paper.
c      Left corner ...... (  0.,  0.)
c      Right corner ..... (600.,780.)
c
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
c
       write(77,*) '90.0 ro'
       write(77,*) '50.0 -550.0 tl'
c
c*  If nfrm=4, four frames in a page (top-left frame).
c
       if(nfrm.eq.1) then
          write(77,*) '1.00 1.00 sc'
       else
          write(77,*) '0.50 0.50 sc'
          write(77,*) '0.0 550.0 tl'
       end if
c
       return
       end
c
c
c-----------------------------
       subroutine gclose
c-----------------------------
       call plote
       return
       end
c
c
c-----------------------------
       subroutine plote
c-----------------------------
       write(77,10) 
   10  format('@eop')
       return
       end
c
c
c-----------------------------------------
       subroutine chart
c-----------------------------------------
c*     Four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
c
c
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
c
c*  Frame 1: open a new page.
c
       if(loc.eq.0) then
          call plote
c
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
c
          write(77,10) 
   10     format('%%PageTrailer    % Need for the page count')
c
          write(77,20) lpage,lpage
   20     format('%%Page:',1x,i2,1x,i2)
c
          write(77,30) 
   30     format('%%BeginPageSetup',/
     *           '%%EndPageSetup',/
     *           '@bop')
c
          write(77,*) '90.0 ro'
          write(77,*) '50.0 -550.0 tl'
c
          if(nfrm.eq.1) then
             write(77,*) '1.00 1.00 sc'
          else
             write(77,*) '0.50 0.50 sc'
             write(77,*) '0.0  550.0 tl'
          end if
c
          return
       end if
c
c
c-----------------------------------------------------
c      First cancel the previous translation, then
c      make a new translation (scale factor alive).
c-----------------------------------------------------
c*   Frames 2-4:
c
       if(loc.eq.1) then
          write(77,*) '  0.0 -550.0 tl'
          write(77,*) '700.0  550.0 tl'
       end if
c
       if(loc.eq.2) then
          write(77,*) '-700.0 -550.0 tl'
          write(77,*) '   0.0    0.0 tl'
       end if
c
       if(loc.eq.3) then
          write(77,*) '  0.0 0.0 tl'
          write(77,*) '700.0 0.0 tl'
       end if
c
       return
       end
c
c
c------------------------------------
       subroutine factor(fct)
c------------------------------------
       write(77,10) fct,fct
   10  format(f6.2,1x,f6.2,' sc')
       return
       end
c
c
c---------------------------------------
       subroutine newcolor (ic,r,g,b)
c---------------------------------------
c  ic= 3 tri-color
c  ic= 0 gray scale, r= 0. for black
c
       write(77,*) 'stroke'
c
       if(ic.eq.0) then
         write(77,10) 1.-r  ! 0. for black
   10    format(f4.1,' setgray')
       end if
c
       if(ic.eq.3) then
         write(77,30) r,g,b
   30    format(3f4.1,' setrgbcolor')
       end if
c
       return
       end
c
c
c------------------------------------
       subroutine newpen (ip)
c------------------------------------
       i1=(ip-1)/2
       i2=ip-2*i1
       write(77,*) 'sn'
       pi1=0.40*float(i1-1)
       write(77,30) pi1
   30  format(f3.1,' sl')
       if(i2.ne.1) then
       write(77,*) '[2 2] 0 sd'
       endif
       return
       end
c
c
c-----------------------------
       subroutine linee
c-----------------------------
       write(77,*) 'st'
       return
       end
c
c
c------------------------------------
       subroutine plot (x0,y0,ip)
c------------------------------------
c
       x= x0
       y= y0
       h= 0.
       n= 777
       call sunscl (x,y,h,n)
c
       if(ip.eq.3)  write(77,10) x,y
       if(ip.eq.2)  write(77,20) x,y
       if(ip.eq.-3) write(77,30) x,y
       if(ip.eq.-2) write(77,40) x,y,x,y
   10  format(f5.1,1x,f5.1,' m')
   20  format(f5.1,1x,f5.1,' l')
   30  format(f5.1,1x,f5.1,' tl')
   40  format(f5.1,1x,f5.1,' l sn',1x,f5.1,1x,f5.1,' tl')
c       write(77,*) 'st'
       return
       end
c
c
c-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
c-------------------------------------------------
       character isymb*80,ica*80,ich(80)*1
       equivalence (ica,ich(1))
c
       x= x0
       y= y0
       h= h0
       n= n0
       call sunscl (x,y,h,n)
c
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
c*
       ica= isymb
       write(77,*) '(',(ich(i),i=1,n),') s'
c
       return
       end
c
c
c-----------------------------------------------
       subroutine number (x0,y0,h0,anu,ang,n0)
c-----------------------------------------------
       character  isymb*9
c
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
c
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
c
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
c
       write(isymb,40) anu
   40  format(1pe9.2)
       write(77,*) '(',isymb,') s'
c
       return
       end
c
c
c---------------------------------------------------
       subroutine sunscl (x,y,h,n)
c---------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
c
       if(x.eq.999.) then
         x= x0 +iabs(n0)*h0
       else
         x= fmag*x
         x0= x
       end if
c
       if(y.eq.999.) then
         y= y0
       else
         y= fmag*y
         y0= y
       end if
c
       h= fmag*h
       h0= h
       if(n.ne.777) n0= n
c
       return
       end
