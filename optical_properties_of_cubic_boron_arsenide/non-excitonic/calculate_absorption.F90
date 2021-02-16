!==========================================================================================================
!PROGRAM calculate_B
!
  !-- Calculate direct optical absorption spectrum, including spin-orbit coupling
  !-- a) Wannier-interpolate eigenvalues and dipole MEs to a fine k-grid

  !-- Required input: u_matrix.dat  u_matrix_opt.dat  prefix.eig  prefix.gw.eig prefix.mmn  b_wb.dat  prefix.wout
  !-- Hard coded parameters (not exhaustive): nw, nb, ivbm, nk, nkint, eta, nomega, omega_min omega_max

  !-- E. Kioupakis, 04/19/2010
  !-- DJB - added inread of Wannier90 output file prefix.wout - 2015.1.21
  !-- DJB - added vme renormalization with GW eigenenergies. Now requires inread of both GW and DFT eigenenergies - 2015.09.02
!===========================================================================================================

PROGRAM main
  IMPLICIT NONE
  INCLUDE 'fftw3.f'
#ifdef MPI
  INCLUDE 'mpif.h'
#endif

  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12,256)
  REAL(r8), PARAMETER :: ryd = 13.605698066_r8, ha2eV = 2.0_r8*ryd, twopi=6.2831853071795864769_r8, bohr2cm =5.29177249E-9_r8, kelvin2au = 3.16682968067E-6_r8, time_au2sec = 2.418884326E-17_r8, bohr2Ang = 0.529177249_r8, pi = 3.14159265358979323846264338327950288_r8


  !-- EDIT to match BAs
  INTEGER, PARAMETER :: nskip = 0 !-- number of bands to skip  
  INTEGER, PARAMETER :: nw    = 30-nskip !-- Number of Wannier functions
  INTEGER, PARAMETER :: nb    = 50-nskip !-- Number of total bands
  INTEGER, PARAMETER :: ivbm  = 4-nskip, icbm = ivbm+1
  INTEGER, PARAMETER :: nk1 = 8, nk2 = 8, nk3 = 8, nktot = nk1*nk2*nk3 !-- Monkhorst-Pack grid
  INTEGER, PARAMETER :: nkint(3) = (/ 128, 128, 128 /) !-- number of k-points to interpolate to, increase for convergence
  INTEGER, PARAMETER :: nomega = 4001 !-- number of frequencies for absorption spectra 
  REAL(r8), PARAMETER :: eta = 0.1/ha2eV !-- broadening of absorption spectra, use a value between 0.01 to 0.1 eV
  REAL(r8), PARAMETER :: omega_min = 0.0_r8/ha2eV, omega_max = 40.0_r8/ha2eV, omega_step = (omega_max-omega_min)/(nomega-1) !-- Range of photon frequencies for absorption




  REAL(r8) :: omega(nomega), absorption_LS(3,nomega), absorption_LS_local(3,nomega), absorption_noLS(3,nomega), absorption_noLS_local(3,nomega)
!  REAL(r8), PARAMETER :: Egap_exp = 6.12101822685397_r8/ha2eV   !-- GW smallest direct gap of h-BN
  REAL(r8) :: Elda(nb,nk1,nk2,nk3), dmin(nk1,nk2,nk3), Egw(nb,nk1,nk2,nk3)
  COMPLEX(r8) :: Umat(nw,nw,nk1,nk2,nk3), Umat_opt(nb,nw,nk1,nk2,nk3), Umat_total(nb,nw,nk1,nk2,nk3)
  INTEGER :: ndimwin(nk1,nk2,nk3)


  !-- With and without spin orbit, L.S, LS
  COMPLEX(r8) :: H0(nb,nb,nk1,nk2,nk3), Hrot(nw,nw,nk1,nk2,nk3), HR(nw,nw,nk1,nk2,nk3), Hint(nw,nw), Htemp(nb,nw) 
  COMPLEX(r8) :: H0gw(nb,nb,nk1,nk2,nk3), Hrotgw(nw,nw,nk1,nk2,nk3), HRgw(nw,nw,nk1,nk2,nk3), Hintgw(nw,nw), Htempgw(nb,nw) 
  COMPLEX(r8) :: H0_LS(nb,nb,2,2,nk1,nk2,nk3), Hrot_LS(nw,nw,2,2,nk1,nk2,nk3), HR_LS(nw,nw,2,2,nk1,nk2,nk3), Hint_undiag(nw,nw,2,2), Hint_LS(2*nw,2*nw)
  COMPLEX(r8) :: H0gw_LS(nb,nb,2,2,nk1,nk2,nk3), Hrotgw_LS(nw,nw,2,2,nk1,nk2,nk3), HRgw_LS(nw,nw,2,2,nk1,nk2,nk3), Hintgw_LS(2*nw,2*nw)
	
  
  !-- Velocity MEs
  COMPLEX(r8) :: Hint_alpha_W(3,nw,nw), Hint_alpha_H(3,nw,nw,nkint(3)), Amn_int_W(3,nw,nw), Amn_int_H(3,nw,nw,nkint(3)),vme(3,nw,nw,nkint(3)),vmegw(3,nw,nw,nkint(3))
  COMPLEX(r8) :: vme_LS(3,2*nw,2*nw,nkint(3)),vmegw_LS(3,2*nw,2*nw,nkint(3)),vme_LS_temp(3,2*nw,2*nw,nkint(3))
  COMPLEX(r8) :: Ame_0(3,nb,nb,nk1,nk2,nk3), Ame_rot(3,nw,nw,nk1,nk2,nk3),Ame_R(3,nw,nw,nk1,nk2,nk3), Ame_int(3,nw,nw,nkint(3)), Ame_int_rot(3,nw,nw,nkint(3)), Ame_temp(3,nb,nw), Mmn_temp(nb,nw)


  REAL(r8) :: output_density, e_cutoff, kvec_min(3), kvec_temp(3), qmin, qtemp
  INTEGER :: nkeep, i_kT
  REAL(r8), ALLOCATABLE :: k_out(:,:)


  REAL(r8), ALLOCATABLE :: Elda_out_temp(:,:), Elda_out_local(:,:), Egw_out_temp(:,:), Egw_out_local(:,:)
  REAL(r8), ALLOCATABLE :: Elda_out_temp_LS(:,:), Elda_out_local_LS(:,:), Egw_out_temp_LS(:,:), Egw_out_local_LS(:,:)

  REAL(r8) :: k_cumm, kdist, k_diff(3), temp_x, temp_y
  COMPLEX(r8), ALLOCATABLE :: Uprime(:,:,:), Uprimegw(:,:,:), Uprime_LS(:,:,:), Uprimegw_LS(:,:,:)

  
  INTEGER :: ndegen(nk1,nk2,nk3) !-- number of equivalent R-vectors in W-S cell for each R vector
  INTEGER :: map(3,125,nk1,nk2,nk3) !-- integer coordinates of each equiv of Rpoint in Wigner-Seitz cell
  REAL(r8) :: dtemp, evbm, f_e, f_h

  REAL(r8) :: Vcell

  REAL(r8) :: d1,d2,d3,d4,d5,d6,a1,a2,a3,a4,a5,a6
  INTEGER :: ik, ikp

  INTEGER :: nbv
  REAL(r8), ALLOCATABLE :: bv(:,:,:), wb(:)
  INTEGER, ALLOCATABLE :: index_bv(:,:,:)
  COMPLEX(r8), ALLOCATABLE :: Mmn(:,:,:), Mmn_rot(:,:,:)


  !-- Interpolate: DOS
  CHARACTER(LEN=10) :: c
  CHARACTER(LEN=100) :: filename, line, keyword, dummy
	
  !-- FFTw:
  COMPLEX(r8) :: fftw_in( NK1,NK2,NK3), fftw_out(NK1,NK2,NK3)
  INTEGER(8) :: plan

  REAL(r8) :: A(3,3), B(3,3), dotprod, tempk(3), tempk_min, tempk_temp(3), tempk_tempmin, kmin(3),tempnorm, tempR(3)
        
  !-- LAPACK
  INTEGER, PARAMETER :: LWORK=4*NW-1
  INTEGER :: INFO
  COMPLEX(r8) :: WORK(LWORK)
  REAL(r8) :: RWORK(14*NW-2), W(2*NW), Wgw(2*NW)

  INTEGER :: i,j,k,k1,k2,k3,m,n,idummy,is,isp,in,inp,isn,isnp,ix,iy,iz,kt1,kt2,kt3,i1,i2,i3,j1,j2,j3, nt, i_loop

  INTEGER nprocs, ident, ierr, myerr, jproc
#ifdef MPI
  INTEGER status(MPI_STATUS_SIZE)

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, ident, ierr)
#else
  nprocs = 1
  ident = 0
#endif


!#ifdef MPI
!  IF (nkint(1)*nkint(2) .NE. nprocs ) THEN
!     PRINT*, 'ERROR: number of processors is equal to nkint(1)*nkint(2)=', nkint(1)*nkint(2)
!     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!     CALL MPI_FINALIZE(ierr)
!  END IF
!#endif


  ALLOCATE(Uprime(nw,nw,nkint(3)))
  ALLOCATE(Uprimegw(nw,nw,nkint(3)))
  ALLOCATE(Uprime_LS(2*nw,2*nw,nkint(3)))
  ALLOCATE(Uprimegw_LS(2*nw,2*nw,nkint(3)))
  ALLOCATE(Elda_out_temp(nw,nkint(3)))
  ALLOCATE(Egw_out_temp(nw,nkint(3)))
  ALLOCATE(Elda_out_local(nw,nkint(3)))
  ALLOCATE(Egw_out_local(nw,nkint(3)))
  ALLOCATE(Elda_out_temp_LS(2*nw,nkint(3)))
  ALLOCATE(Egw_out_temp_LS(2*nw,nkint(3)))
  ALLOCATE(Elda_out_local_LS(2*nw,nkint(3)))
  ALLOCATE(Egw_out_local_LS(2*nw,nkint(3)))


  Elda_out_local = 0.0_r8
  Egw_out_local = 0.0_r8

  vme = CMPLX(0.0_r8,0.0_r8)
  vmegw = CMPLX(0.0_r8,0.0_r8)
  vme_LS = CMPLX(0.0_r8,0.0_r8)
  vmegw_LS = CMPLX(0.0_r8,0.0_r8)
  
  !-- EDIT read parameters from Wannier output file, specify filename
  !-- Read in Volume, Real & Reciprocal lattice vectors
  IF (ident .EQ. 0 ) PRINT*, 'Reading file: prefix.wout'
  OPEN( UNIT = 101, FILE = 'BAs.wout', ACTION = 'READ' )



  DO WHILE( .TRUE. )
     READ( 101, '(a100)' ) line
     line = ADJUSTL( line )
     keyword = line( 1 : SCAN( line, ' ' ) - 1 )
     line = ADJUSTL( line( SCAN( line, ' ' ) + 1 : 100 ) )
     IF ( TRIM(keyword) .EQ. 'Lattice' ) THEN
        DO j = 1, 3
           READ( 101, * ) dummy, A( 1, j ), A( 2, j ), A( 3, j )
        END DO
        A = A / bohr2Ang
     ELSEIF ( TRIM(keyword) .EQ. 'Unit' ) THEN
        READ( line, * ) dummy, dummy, Vcell
        Vcell = Vcell / ( bohr2Ang * bohr2Ang * bohr2Ang )
     ELSEIF ( TRIM(keyword) .EQ. 'Reciprocal-Space' ) THEN
        DO j = 1, 3
           READ( 101, * ) dummy, B( 1, j ), B( 2, j ), B( 3, j )
        END DO
        B = B * bohr2Ang
        EXIT
     END IF
  END DO
  CLOSE(101)
	
  !-- First pass: find equivalent R-point in 1st BZ and minimum distance from origin
  IF ( ident .EQ. 0 ) PRINT*, 'First pass: equivalent R-points'
  DO k1 = 1, NK1; DO k2 = 1, NK2; DO k3 = 1, NK3
     dmin(k1,k2,k3) = SQRT( ((k1-1)*A(1,1)+(k2-1)*A(1,2)+(k3-1)*A(1,3))**2 + &
                          & ((k1-1)*A(2,1)+(k2-1)*A(2,2)+(k3-1)*A(2,3))**2 + &
                          & ((k1-1)*A(3,1)+(k2-1)*A(3,2)+(k3-1)*A(3,3))**2)
     map(:,1,k1,k2,k3) = (/ k1-1, k2-1, k3-1 /)
           
     DO ix = -2, 2; DO iy = -2, 2; DO iz = -2, 2
        kt1 = k1 + ix * NK1 -1 
        kt2 = k2 + iy * NK2 -1
        kt3 = k3 + iz * NK3 -1
        dtemp = SQRT( (kt1*A(1,1)+kt2*A(1,2)+kt3*A(1,3))**2 + &
                    & (kt1*A(2,1)+kt2*A(2,2)+kt3*A(2,3))**2 + &
                    & (kt1*A(3,1)+kt2*A(3,2)+kt3*A(3,3))**2 )
        IF ( dtemp .LT. dmin(k1,k2,k3) ) THEN
           dmin(k1,k2,k3) = dtemp
           map(:,1,k1,k2,k3) = (/ kt1, kt2, kt3 /)
        END IF
     END DO; END DO; END DO
     !PRINT*, k1,k2,k3,map(1,k1,k2,k3),map(2,k1,k2,k3),map(3,k1,k2,k3),dmin
	   
  END DO; END DO; END DO

  !-- Second pass: find degeneracies:
  IF (ident .EQ. 0 ) PRINT*, 'Second pass: equivalent R-points'
  DO k1 = 1, NK1
  DO k2 = 1, NK2
  DO k3 = 1, NK3
     ndegen(k1,k2,k3) = 0 !-- 0 because we should find the original point again
     DO ix = -2, 2
     DO iy = -2, 2
     DO iz = -2, 2
        kt1 = k1 + ix * NK1 -1 
        kt2 = k2 + iy * NK2 -1
        kt3 = k3 + iz * NK3 -1
        dtemp = SQRT( (kt1*A(1,1)+kt2*A(1,2)+kt3*A(1,3))**2 + &
                    & (kt1*A(2,1)+kt2*A(2,2)+kt3*A(2,3))**2 + &
                    & (kt1*A(3,1)+kt2*A(3,2)+kt3*A(3,3))**2 )
        IF ( ABS(dtemp-dmin(k1,k2,k3)) .LT. 1.0E-4 ) THEN !-- the two points are equivalent
           ndegen(k1,k2,k3) = ndegen(k1,k2,k3) + 1
           map(1,ndegen(k1,k2,k3),k1,k2,k3) = kt1
           map(2,ndegen(k1,k2,k3),k1,k2,k3) = kt2
           map(3,ndegen(k1,k2,k3),k1,k2,k3) = kt3
        END IF
     END DO
     END DO
     END DO
     !DO i = 1, ndegen(k1,k2,k3)
     !   WRITE(*,("7i4,F12.8")) k1,k2,k3,ndegen(k1,k2,k3),map(:,i,k1,k2,k3),dmin(k1,k2,k3)
     !END DO
  END DO
  END DO
  END DO

  IF (ident .EQ. 0 ) PRINT*, 'Reading files: b_wb.dat'
  OPEN(UNIT=1,FILE='b_wb.dat',ACTION='READ')
  READ(1,*) nbv
  ALLOCATE(bv(3,nbv,nktot),index_bv(3,nbv,nktot),wb(nbv),Mmn(nb,nb,nbv), Mmn_rot(nw,nw,nbv))
  DO i = 1, nktot
     DO j = 1, nbv
        READ(1,*) k1, k2, bv(:,j,i)
     END DO
  END DO
  DO j = 1, nbv
     READ(1,*) k1, wb(j)
  END DO
  !write(*,*) wb
  CLOSE(1)

  !-- Convert from Ang-1 to Bohr units:
  bv = bv * bohr2Ang
  wb = wb / bohr2Ang**2


  IF (ident .EQ. 0 ) PRINT*, 'Reading files: u_matrix_opt.dat'
  OPEN(UNIT=1,FILE='u_matrix_opt.dat',ACTION='READ')
  READ(1,*)
  Umat_opt = CMPLX(0.0_r8,0.0_r8)
  DO k1 = 1, NK1
  DO k2 = 1, NK2
  DO k3 = 1, NK3
     READ(1,*) ndimwin(k1,k2,k3)
     DO m = 1, NW
     DO n = 1, ndimwin(k1,k2,k3)
        READ(1,*) idummy,idummy,idummy,temp_x,temp_y
        Umat_opt(n,m,k1,k2,k3) = CMPLX(temp_x,temp_y)
     END DO
     END DO
  END DO
  END DO
  END DO
  CLOSE(1)

  IF (ident .EQ. 0 )PRINT*, 'Reading files: u_matrix.dat'
  OPEN(UNIT=1,FILE='u_matrix.dat',ACTION='READ')
  READ(1,*)
  Umat = CMPLX(0.0_r8,0.0_r8)
  DO k1 = 1, NK1
  DO k2 = 1, NK2
  DO k3 = 1, NK3
     DO m = 1, NW
     DO n = 1, NW
        READ(1,*) idummy,idummy,idummy,temp_x,temp_y
        Umat(n,m,k1,k2,k3) = CMPLX(temp_x,temp_y)
     END DO
     END DO
  END DO
  END DO
  END DO
  CLOSE(1)


  !-- EDIT specify if reading LDA or GW eigenvalues
  IF (ident .EQ. 0 )PRINT*, 'Reading files: .eig'
  OPEN(UNIT=1,FILE='BAs.LDA.eig',ACTION='READ')
  DO k1 = 1, NK1
  DO k2 = 1, NK2
  DO k3 = 1, NK3
     DO in = 1, NB
        READ(1,*) idummy, idummy, elda(in,k1,k2,k3)
        !-- convert to Hartree
        elda(in,k1,k2,k3) = elda(in,k1,k2,k3)/ha2eV
     END DO
  END DO
  END DO
  END DO
  CLOSE(1)

  IF (ident .EQ. 0 )PRINT*, 'Reading files: .gw.eig'
  OPEN(UNIT=1,FILE='BAs.GW.eig',ACTION='READ')
  DO k1 = 1, NK1
  DO k2 = 1, NK2
  DO k3 = 1, NK3
     DO in = 1, NB
        READ(1,*) idummy, idummy, egw(in,k1,k2,k3)
        !-- convert to Hartree
        egw(in,k1,k2,k3) = egw(in,k1,k2,k3)/ha2eV
     END DO
  END DO
  END DO
  END DO
  CLOSE(1)
  IF (ident .EQ. 0 )PRINT*, 'Done.'


  !-- EDIT: use the eigenvalues for which you ran wannier. I think you used the GW eigenvalues when running wannier90 for BAs, right? then keep this line uncommented
  elda = egw


  !!from Kyle - editing to not use somatrix
  IF (ident .EQ. 0 )PRINT*, 'Reading files: somatrix'
  !OPEN(unit=1,file='somatrix.dat',ACTION='READ') 
  DO k1 = 1, NK1;  DO k2 = 1, NK2;  DO k3 = 1, NK3
     !!READ(1,*)
     !!READ(1,*)
     DO is  = 1, 2;     DO isp = 1, 2
        DO in  = 1, NB;        DO inp = 1, NB
           !!read(1,*) idummy,idummy,idummy,idummy,temp_x,temp_y


           !!! TURN OFF SPIN ORBIT< DEBUG
           !!H0_LS(in,inp,is,isp,k1,k2,k3) = CMPLX(temp_x,temp_y)*ryd/ha2eV
           H0_LS(in,inp,is,isp,k1,k2,k3) = CMPLX(0.0_r8,0.0_r8)*ryd/ha2eV !!!--- EDIT uncomment this line to turn off spin-orbit


           H0gw_LS(in,inp,is,isp,k1,k2,k3) = H0_LS(in,inp,is,isp,k1,k2,k3)
           !             H0(in,inp,is,isp,k1,k2,k3) = CMPLX(0.0D0,0.0D0)  
           IF ( is .EQ. isp .AND. in .EQ. inp ) THEN

!!!--- EDIT: if wannierization was done with GW eigenvalues keep the next lines unchanged
              H0_LS(in,inp,is,isp,k1,k2,k3) = H0_LS(in,inp,is,isp,k1,k2,k3) + Elda(in,k1,k2,k3) !-- Wannier done with GW
              !H0_LS(in,inp,is,isp,k1,k2,k3) = H0_LS(in,inp,is,isp,k1,k2,k3) + Egw(in,k1,k2,k3)   !-- wannier done with LDA




              H0gw_LS(in,inp,is,isp,k1,k2,k3) = H0gw_LS(in,inp,is,isp,k1,k2,k3) + Egw(in,k1,k2,k3)
           END IF
        END DO;        END DO
     END DO;     END DO
  END DO;  END DO;  END DO
  !!CLOSE(1)


  H0 = CMPLX(0.0_r8,0.0_r8)
  H0gw = CMPLX(0.0_r8,0.0_r8)
  DO k1 = 1, NK1;  DO k2 = 1, NK2;  DO k3 = 1, NK3
     DO in  = 1, NB;        DO inp = 1, NB
        IF ( in .EQ. inp ) THEN
           H0(in,inp,k1,k2,k3) = CMPLX(Elda(in,k1,k2,k3),0.0_r8)
           H0gw(in,inp,k1,k2,k3) = CMPLX(Egw(in,k1,k2,k3),0.0_r8)
        END IF
     END DO;        END DO
  END DO;  END DO;  END DO



  !-- rotate:
  IF (ident .EQ. 0 ) PRINT*, 'Rotating...'

  !-- Combine Umat-rices
  DO k1 = 1, NK1
  DO k2 = 1, NK2
  DO k3 = 1, NK3
     Umat_total(:,:,k1,k2,k3) = MATMUL(Umat_opt(:,:,k1,k2,k3),Umat(:,:,k1,k2,k3))
  END DO
  END DO 
  END DO

  IF (ident .EQ. 0 ) THEN
     PRINT*, 'Reading file MMN...'

     OPEN(UNIT=1,FILE='BAs.mmn',ACTION='READ')
     READ(1,*)
     READ(1,*)
     Ame_rot = CMPLX(0.0_r8,0.0_r8)
     DO i = 1, nktot
        DO k = 1, nbv
           READ(1,*) ik, ikp, i1,i2,i3
           k1 = (ik-1)/(nk2*nk3) + 1
           k2 = (ik-1-(k1-1)*(nk2*nk3))/nk3 + 1
           k3 = MOD(ik-1,nk3)+1
           !PRINT*, i, k1,k2,k3
           kt1 = (ikp-1)/(nk2*nk3) + 1
           kt2 = (ikp-1-(kt1-1)*(nk2*nk3))/nk3 + 1
           kt3 = MOD(ikp-1,nk3)+1
           !WRITE(*,("8I5")) ik, k1,k2,k3,ikp, kt1,kt2,kt3

           DO m = 1, nb
           DO n = 1, nb
              READ(1,*) a1,a2
              Mmn(n,m,k) = CMPLX(a1,a2)
           END DO
           END DO
           
           !-- Rotate Mmn:
           Mmn_temp(1:ndimwin(k1,k2,k3),1:NW) = MATMUL( Mmn(1:ndimwin(k1,k2,k3),1:ndimwin(kt1,kt2,kt3),k), & 
                & Umat_total(1:ndimwin(kt1,kt2,kt3),1:NW,kt1,kt2,kt3))
           Mmn_rot(1:NW,1:NW,k) = MATMUL( TRANSPOSE(CONJG(Umat_total(1:ndimwin(k1,k2,k3),1:NW,k1,k2,k3))),Mmn_temp(1:ndimwin(k1,k2,k3),1:NW))
           
           !-- subtract delta_mn
           DO m = 1, NW
              Mmn_rot(m,m,k) = Mmn_rot(m,m,k) - 1.0_r8
           END DO
           DO j = 1, 3
              Ame_rot(j,:,:,k1,k2,k3) = Ame_rot(j,:,:,k1,k2,k3) + CMPLX(0.0_r8,1.0_r8) * wb(k)*bv(j,k,ik)*Mmn_rot(:,:,k)
           END DO
        END DO
     END DO
     CLOSE(1)
     PRINT*, 'DONE: Reading file MMN.'
  END IF

#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST( Ame_rot, 3*nw*nw*nk1*nk2*nk3, MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif


  IF (ident .EQ. 0 ) PRINT*, 'Rotating eigenvalues, pmes...'
  DO k1 = 1, NK1; DO k2 = 1, NK2; DO k3 = 1, NK3

     Htemp(1:ndimwin(k1,k2,k3),1:nw) = MATMUL(H0(1:ndimwin(k1,k2,k3),1:ndimwin(k1,k2,k3),k1,k2,k3),& 
          & Umat_total(1:ndimwin(k1,k2,k3),1:nw,k1,k2,k3))
     Hrot(1:nw,1:nw,k1,k2,k3)  = MATMUL( TRANSPOSE(CONJG(Umat_total(1:ndimwin(k1,k2,k3),1:nw,k1,k2,k3))), & 
          & Htemp(1:ndimwin(k1,k2,k3),1:nw))

     Htempgw(1:ndimwin(k1,k2,k3),1:nw) = MATMUL(H0gw(1:ndimwin(k1,k2,k3),1:ndimwin(k1,k2,k3),k1,k2,k3),& 
          & Umat_total(1:ndimwin(k1,k2,k3),1:nw,k1,k2,k3))
     Hrotgw(1:nw,1:nw,k1,k2,k3)  = MATMUL( TRANSPOSE(CONJG(Umat_total(1:ndimwin(k1,k2,k3),1:nw,k1,k2,k3))), & 
          & Htempgw(1:ndimwin(k1,k2,k3),1:nw))

     DO is  = 1, 2;    DO isp = 1, 2
     Htemp(1:ndimwin(k1,k2,k3),1:nw) = MATMUL(H0_LS(1:ndimwin(k1,k2,k3),1:ndimwin(k1,k2,k3),is,isp,k1,k2,k3),& 
          & Umat_total(1:ndimwin(k1,k2,k3),1:nw,k1,k2,k3))
     Hrot_LS(1:nw,1:nw,is,isp,k1,k2,k3)  = MATMUL( TRANSPOSE(CONJG(Umat_total(1:ndimwin(k1,k2,k3),1:nw,k1,k2,k3))), & 
          & Htemp(1:ndimwin(k1,k2,k3),1:nw))

     Htempgw(1:ndimwin(k1,k2,k3),1:nw) = MATMUL(H0gw_LS(1:ndimwin(k1,k2,k3),1:ndimwin(k1,k2,k3),is,isp,k1,k2,k3),& 
          & Umat_total(1:ndimwin(k1,k2,k3),1:nw,k1,k2,k3))
     Hrotgw_LS(1:nw,1:nw,is,isp,k1,k2,k3)  = MATMUL( TRANSPOSE(CONJG(Umat_total(1:ndimwin(k1,k2,k3),1:nw,k1,k2,k3))), & 
          & Htempgw(1:ndimwin(k1,k2,k3),1:nw))
     END DO; END DO
  END DO; END DO; END DO
  IF (ident .EQ. 0 ) PRINT*, 'Done rotating!'

  !-- FFT
  IF (ident .EQ. 0 ) PRINT*, 'FFTing...'
  CALL dfftw_plan_dft_3d(plan,NK1,NK2,NK3,fftw_in,fftw_out,FFTW_FORWARD,FFTW_MEASURE)
  DO in  = 1, NW
  DO inp = 1, NW
     fftw_in(:,:,:) = Hrot(in,inp,:,:,:)
     CALL dfftw_execute(plan)
     HR(in,inp,:,:,:) = fftw_out(:,:,:)
  END DO
  END DO
  DO is  = 1, 2;  DO isp = 1, 2
  DO in  = 1, NW
  DO inp = 1, NW
     fftw_in(:,:,:) = Hrot_LS(in,inp,is,isp,:,:,:)
     CALL dfftw_execute(plan)
     HR_LS(in,inp,is,isp,:,:,:) = fftw_out(:,:,:)
  END DO
  END DO
  END DO; END DO

  IF (ident .EQ. 0 ) PRINT*, 'FFTing GW...'
  DO in  = 1, NW
  DO inp = 1, NW
     fftw_in(:,:,:) = Hrotgw(in,inp,:,:,:)
     CALL dfftw_execute(plan)
     HRgw(in,inp,:,:,:) = fftw_out(:,:,:)
  END DO
  END DO
  DO is  = 1, 2;  DO isp = 1, 2
  DO in  = 1, NW
  DO inp = 1, NW
     fftw_in(:,:,:) = Hrotgw_LS(in,inp,is,isp,:,:,:)
     CALL dfftw_execute(plan)
     HRgw_LS(in,inp,is,isp,:,:,:) = fftw_out(:,:,:)
  END DO
  END DO
  END DO; END DO

  DO in  = 1, NW
  DO inp = 1, NW
     DO j = 1, 3
        fftw_in(:,:,:) = ame_rot(j,in,inp,:,:,:)
        CALL dfftw_execute(plan)
        ame_R(j,in,inp,:,:,:) = fftw_out(:,:,:)
     END DO
  END DO
  END DO
  IF (ident .EQ. 0 ) PRINT*, 'Done FFTing!'

  IF ( ident .EQ. 0 ) THEN
  OPEN(UNIT=1,FILE='HmnR_vs_R.dat')
  DO in  = 1, NW
  DO inp = 1, NW
     DO k1 = 1, NK1
     DO k2 = 1, NK2
     DO k3 = 1, NK3
        tempR(:) = map(1,1,k1,k2,k3)*A(:,1)+ map(2,1,k1,k2,k3)*A(:,2)+ map(3,1,k1,k2,k3)*A(:,3)
        WRITE(1,"(2F20.8)") SQRT(DOT_PRODUCT(tempR,tempR)), REAL(HR(in,inp,k1,k2,k3)*CONJG(HR(in,inp,k1,k2,k3)))
     END DO
     END DO
     END DO
  END DO
  END DO
  CLOSE(1)


  OPEN(UNIT=1,FILE='Ame_vs_R.dat')
  DO in  = 1, NW
  DO inp = 1, NW
     DO k1 = 1, NK1
     DO k2 = 1, NK2
     DO k3 = 1, NK3
        tempR(:) = map(1,1,k1,k2,k3)*A(:,1)+ map(2,1,k1,k2,k3)*A(:,2)+ map(3,1,k1,k2,k3)*A(:,3)
        WRITE(1,"(4F20.8)") SQRT(DOT_PRODUCT(tempR,tempR)), REAL(Ame_R(:,in,inp,k1,k2,k3)*CONJG(Ame_R(:,in,inp,k1,k2,k3)))
     END DO
     END DO
     END DO
  END DO
  END DO
  CLOSE(1)

  END IF

  IF (ident .EQ. 0 ) PRINT*, 'Interpolating, direct...'

!-- Approach 1:
!  !DO i1 = 1, nkint(1) 
!  DO i1 = ident+1, nkint(1), nprocs
!     !WRITE(*,"(I3,A5,I3)")  i1, ' of ', nkint(1)
!  DO i2 = 1, nkint(2) 
!     WRITE(*,*) 'i2=', i2


  DO i = 1, nomega
     omega(i) = omega_min + (i-1)*omega_step
  END DO
  absorption_LS = 0.0_r8
  absorption_LS_local = 0.0_r8
  absorption_noLS = 0.0_r8
  absorption_noLS_local = 0.0_r8

!-- Approach 2:
  iloop: DO i_loop = ident, nkint(1)*nkint(2)-1, nprocs
     i1 = i_loop/nkint(2) + 1
     i2 = mod(i_loop,nkint(2)) + 1
  DO i3 = 1, nkint(3)


!  k_extrema(:,1) = (/ 0.0_r8, 0.0_r8, 0.0_r8 /) !-- direct gap
!  k_extrema(:,2) = (/ 0.495370_r8, 0.5_r8, 0.185185_r8 /) !-- CBM
!  k_extrema(:,3) = (/ 0.5_r8, 0.5_r8, 0.055556_r8 /)  !-- min direct
!  DO i3 = 1, 1
     IF (ident .EQ. 0 ) WRITE(*,*) 'i3=', i3
     !WRITE(*,"(3I3,A5,3I3)")  i1, i2, i3, ' of ', nkint(1), nkint(2), nkint(3)

     tempk = (/ REAL(i1-1,r8)/nkint(1), REAL(i2-1,r8)/nkint(2), REAL(i3-1,r8)/nkint(3) /)
     !tempk(:) = k_extrema(:,i3)

     !-- Find equivalent interpolated kpoint closest to the origin
     tempk_min = SQRT(DOT_PRODUCT(MATMUL(B,tempk), MATMUL(B,tempk)))
     DO ix = -2, 2; DO iy = -2, 2; DO iz = -2, 2
        tempk_temp = (/ REAL(i1-1+ix*nkint(1),r8)/nkint(1),REAL(i2-1+iy*nkint(2),r8)/nkint(2),REAL(i3-1+iz*nkint(3),r8)/nkint(3) /)
        tempk_tempmin = SQRT(DOT_PRODUCT(MATMUL(B,tempk_temp), MATMUL(B,tempk_temp)))
        IF ( tempk_tempmin .LT. tempk_min ) THEN
           tempk_min = tempk_tempmin
          tempk = tempk_temp
        END IF
     END DO; END DO; END DO

     !WRITE(*,("I4,6F10.5")) i3, (/ REAL(i1-1,r8)/nkint(1), REAL(i2-1,r8)/nkint(2), REAL(i3-1,r8)/nkint(3) /), tempk

     Hint = CMPLX(0.0_r8,0.0_r8)
     Hintgw = CMPLX(0.0_r8,0.0_r8)

     DO in  = 1, NW
     DO inp = 1, NW
        DO k1 = 1, NK1
        DO k2 = 1, NK2
        DO k3 = 1, NK3
           DO i = 1, ndegen(k1,k2,k3)
              dotprod = twopi* DOT_PRODUCT(tempk,map(:,i,k1,k2,k3))
              temp_x = COS(dotprod)
              temp_y = SIN(dotprod) !sign????
              Hint(in,inp) = Hint(in,inp) + HR(in,inp,k1,k2,k3) * CMPLX(temp_x,temp_y)/(NKTOT*ndegen(k1,k2,k3))
              Hintgw(in,inp) = Hintgw(in,inp) + HRgw(in,inp,k1,k2,k3) * CMPLX(temp_x,temp_y)/(NKTOT*ndegen(k1,k2,k3))
           END DO
        END DO
        END DO
        END DO
     END DO
     END DO
     

     Hint_LS = CMPLX(0.0_r8,0.0_r8)
     Hintgw_LS = CMPLX(0.0_r8,0.0_r8)
     
     Hint_undiag = CMPLX(0.0_r8,0.0_r8)

     !DO is = 1, 2
     !DO isp = 1, 2
     !DO in  = 1, NW
     !DO inp = 1, NW
        !isn  = in  + (is -1)*NW
        !isnp = inp + (isp-1)*NW

        DO k1 = 1, NK1
        DO k2 = 1, NK2
        DO k3 = 1, NK3
           DO i = 1, ndegen(k1,k2,k3)
              dotprod = twopi* DOT_PRODUCT(tempk,map(:,i,k1,k2,k3))
              temp_x = COS(dotprod)
              temp_y = SIN(dotprod) !sign????
              !Hint_LS(isn,isnp) = Hint_LS(isn,isnp) + HR_LS(in,inp,is,isp,k1,k2,k3) * CMPLX(temp_x,temp_y)/(NKTOT*ndegen(k1,k2,k3))
              !Hintgw_LS(isn,isnp) = Hintgw_LS(isn,isnp) + HRgw_LS(in,inp,is,isp,k1,k2,k3) * CMPLX(temp_x,temp_y)/(NKTOT*ndegen(k1,k2,k3))
              Hint_undiag(:,:,:,:) = Hint_undiag(:,:,:,:) + HR_LS(:,:,:,:,k1,k2,k3) * CMPLX(temp_x,temp_y)/(NKTOT*ndegen(k1,k2,k3))

           END DO
        END DO
        END DO
        END DO
     !END DO
     !END DO
     !END DO
     !END DO
	   
     !--diagonalize:
     !-- save eigenvectors, to be used to rotate interpolated matrices later
     CALL zheev('V','U',NW,Hint,NW,W,WORK,LWORK,RWORK,INFO)
     IF ( INFO .NE. 0 ) THEN
        PRINT*, 'ERROR: diagonalizing'
        STOP
     END IF
     Elda_out_local(1:nw,i3) = W(1:nw)
     Uprime(1:nw,1:nw,i3) = Hint(1:nw,1:nw)
     !PRINT*, Elda_out_local(1:nw,i3)*ha2eV


     !-- rotate spin-orbit hamiltonian on the same basis (same as velocity MEs

     DO is = 1, 2
     DO isp = 1, 2
        Hint_undiag(:,:,is,isp) = MATMUL( Hint_undiag(:,:,is,isp), Uprime(:,:,i3))
        Hint_undiag(:,:,is,isp) = MATMUL( TRANSPOSE(CONJG(Uprime(:,:,i3))),Hint_undiag(:,:,is,isp))
     END DO
     END DO

     DO is = 1, 2
     DO isp = 1, 2
     DO in  = 1, NW
     DO inp = 1, NW
        isn  = in  + (is -1)*NW
        isnp = inp + (isp-1)*NW

        Hint_LS(isn,isnp) = Hint_undiag(in,inp,is,isp)

     END DO
     END DO
     END DO
     END DO

     CALL zheev('V','U',2*NW,Hint_LS,2*NW,W,WORK,LWORK,RWORK,INFO)
     IF ( INFO .NE. 0 ) THEN
        PRINT*, 'ERROR: diagonalizing'
        STOP
     END IF
     Elda_out_local_LS(:,i3) = W(:)
     Uprime_LS(:,:,i3) = Hint_LS(:,:)
     !PRINT*, Elda_out_local_LS(1:2*nw,i3)*ha2eV

     !STOP
     CALL zheev('V','U',NW,Hintgw,NW,Wgw,WORK,LWORK,RWORK,INFO)
     Egw_out_local(1:nw,i3) = Wgw(1:nw)
     Uprimegw(1:nw,1:nw,i3) = Hintgw(1:nw,1:nw)

     CALL zheev('V','U',2*NW,Hintgw_LS,2*NW,Wgw,WORK,LWORK,RWORK,INFO)
     Egw_out_local_LS(:,i3) = Wgw(:)
     Uprimegw_LS(:,:,i3) = Hintgw_LS(:,:)



     Hint_alpha_W = CMPLX(0.0_r8,0.0_r8)
     DO in  = 1, NW
     DO inp = 1, NW
        DO k1 = 1, NK1
        DO k2 = 1, NK2
        DO k3 = 1, NK3
           DO i = 1, ndegen(k1,k2,k3)
              dotprod = twopi* DOT_PRODUCT(tempk,map(:,i,k1,k2,k3))
              temp_x = COS(dotprod)
              temp_y = SIN(dotprod) !sign????
              Hint_alpha_W(:,in,inp) = Hint_alpha_W(:,in,inp) + CMPLX(0.0_r8,1.0_r8)*MATMUL(A,map(:,i,k1,k2,k3))* &
                   & HR(in,inp,k1,k2,k3) * &
                   & CMPLX(temp_x,temp_y)/(NKTOT*ndegen(k1,k2,k3))
           END DO
        END DO
        END DO
        END DO
     END DO
     END DO

     !-- Use LDA eigenvectors
     DO j = 1, 3
        Hint_alpha_H(j,:,:,i3) = MATMUL( Hint_alpha_W(j,:,:), Uprime(:,:,i3))
        Hint_alpha_H(j,:,:,i3) = MATMUL( TRANSPOSE(CONJG(Uprime(:,:,i3))),Hint_alpha_H(j,:,:,i3))
     END DO

     Amn_int_W = CMPLX(0.0_r8,0.0_r8)
     DO in  = 1, NW
     DO inp = 1, NW
        DO k1 = 1, NK1
        DO k2 = 1, NK2
        DO k3 = 1, NK3
           DO i = 1, ndegen(k1,k2,k3)
              dotprod = twopi* DOT_PRODUCT(tempk,map(:,i,k1,k2,k3))
              temp_x = COS(dotprod)
              temp_y = SIN(dotprod) !sign????
              Amn_int_W(:,in,inp) = Amn_int_W(:,in,inp) + &
                   & Ame_R(:,in,inp,k1,k2,k3) * CMPLX(temp_x,temp_y)/(NKTOT*ndegen(k1,k2,k3))
           END DO
        END DO
        END DO
        END DO
     END DO
     END DO

     DO j = 1, 3
        Amn_int_H(j,:,:,i3) = MATMUL( Amn_int_W(j,:,:), Uprime(:,:,i3))
        Amn_int_H(j,:,:,i3) = MATMUL( TRANSPOSE(CONJG(Uprime(:,:,i3))),Amn_int_H(j,:,:,i3))
     END DO

     DO n = 1, nw
     DO m = 1, nw
        vme(:,m,n,i3) = Hint_alpha_H(:,m,n,i3) -CMPLX(0.0_r8,1.0_r8)*(Elda_out_local(n,i3)-Elda_out_local(m,i3))*Amn_int_H(:,m,n,i3)  !-- DJB: this is not formally consistent with our renormalization method, but we're not changing it for now, correct?
        vmegw(:,m,n,i3) = vme(:,m,n,i3)
     END DO
     END DO

     DO n = icbm, nw
     DO m = 1, ivbm
        vmegw(:,m,n,i3) = ((Egw_out_local(n,i3)-Egw_out_local(m,i3))/(Elda_out_local(n,i3)-Elda_out_local(m,i3))) * &
             & vmegw(:,m,n,i3)
        vmegw(:,n,m,i3) = ((Egw_out_local(n,i3)-Egw_out_local(m,i3))/(Elda_out_local(n,i3)-Elda_out_local(m,i3))) * &
             & vmegw(:,n,m,i3)
     END DO
     END DO


     vme_LS(:,:,:,i3) = CMPLX(0.0_r8,0.0_r8)
     vme_LS_temp(:,:,:,i3) = CMPLX(0.0_r8,0.0_r8)
     vmegw_LS(:,:,:,i3) = CMPLX(0.0_r8,0.0_r8)


!     !-- Undo rotation?
!     DO j = 1, 3
!        vme(j,:,:,i3) = MATMUL(vme(j,:,:,i3),TRANSPOSE(CONJG(Uprime(:,:,i3))))
!        vme(j,:,:,i3) = MATMUL(Uprime(:,:,i3),vme(j,:,:,i3))
!     END DO

     DO is = 1, 2
     DO isp = 1, 2
     DO in  = 1, NW
     DO inp = 1, NW
        isn  = in  + (is -1)*NW
        isnp = inp + (isp-1)*NW
        
        !isn = is + (in-1)*2
        !isnp = isp + (inp-1)*2

        IF (is .EQ. isp ) THEN

           vme_LS_temp(:,isn,isnp,i3) = vme(:,in,inp,i3)
        END IF


     END DO
     END DO
     END DO
     END DO 

!     vme_LS(1,:,:,i3) = MATMUL( TRANSPOSE(CONJG(Uprime_LS(:,:,i3))), Uprime_LS(:,:,i3) )
!     DO i = 1, 2*NW
!     DO j = 1, 2*NW
!        PRINT*, i, j, vme_LS(1,i,j,i3)
!     END DO
!     END DO   
      



  
     DO j = 1, 3
        vme_LS(j,:,:,i3) = MATMUL( vme_LS_temp(j,:,:,i3), Uprime_LS(:,:,i3))
        vme_LS(j,:,:,i3) = MATMUL( TRANSPOSE(CONJG(Uprime_LS(:,:,i3))), vme_LS(j,:,:,i3))


        !!!! DEBUG
        !vme_LS(j,:,:,i3) = MATMUL( TRANSPOSE(CONJG(Uprime_LS(:,:,i3))), MATMUL(vme_LS_temp(j,:,:,i3), Uprime_LS(:,:,i3)))

     END DO

     vmegw_LS(:,:,:,i3) = vme_LS(:,:,:,i3)
     DO n = 2*ivbm+1, 2*nw
     DO m = 1, 2*ivbm
        vmegw_LS(:,m,n,i3) = ((Egw_out_local_LS(n,i3)-Egw_out_local_LS(m,i3))/(Elda_out_local_LS(n,i3)-Elda_out_local_LS(m,i3))) * &
             & vmegw_LS(:,m,n,i3)
        vmegw_LS(:,n,m,i3) = ((Egw_out_local_LS(n,i3)-Egw_out_local_LS(m,i3))/(Elda_out_local_LS(n,i3)-Elda_out_local_LS(m,i3))) * &
             & vmegw_LS(:,n,m,i3)
     END DO
     END DO

     !PRINT*, 'No LS: '
     !DO i = ivbm-5, ivbm
     !DO j = ivbm+1, ivbm+6
     ! !WRITE(*,"(2I3,7E20.10)") i,j, vmegw(:,i,j,i3), REAL(DOT_PRODUCT(vmegw(:,i,j,i3),vmegw(:,i,j,i3)))
     !   WRITE(*,"(2I3,2E20.10)") i,j, REAL(DOT_PRODUCT(vme(:,i,j,i3),vme(:,i,j,i3))), REAL(DOT_PRODUCT(vmegw(:,i,j,i3),vmegw(:,i,j,i3)))
     !END DO
     !END DO

     !PRINT*, 'With LS: '
     !DO i = 2*ivbm-11, 2*ivbm
     !DO j = 2*ivbm+1, 2*ivbm+12
     !   !WRITE(*,"(2I3,7E20.10)") i,j, vmegw_LS(:,i,j,i3), REAL(DOT_PRODUCT(vmegw_LS(:,i,j,i3),vmegw_LS(:,i,j,i3)))
     !   WRITE(*,"(2I3,2E20.10)") i,j, REAL(DOT_PRODUCT(vme_LS(:,i,j,i3),vme_LS(:,i,j,i3))), REAL(DOT_PRODUCT(vmegw_LS(:,i,j,i3),vmegw_LS(:,i,j,i3)))
     !END DO
     !END DO



     !IF (i3 .EQ. 1 ) THEN
     !   PRINT*, 'VBM'
     !ELSE IF (i3 .EQ. 2 ) THEN
     !   PRINT*, 'CBM'
     !ELSE 
     !   PRINT*, 'min direct'
     !END IF
     !PRINT*, EGW_out_local(ivbm,i3)*ha2eV, EGW_out_local(ivbm+1,i3)*ha2eV
     !!PRINT*, Elda_out_local(:,i3)*ha2eV
     !PRINT*, SQRT(DOT_PRODUCT(vmegw(:,8,9,i3),vmegw(:,8,9,i3))),  vmegw(:,8,9,i3)


     !-- Calculate absorption:



     !-- EDIT Debug: print matrix elements
     !PRINT*, REAL(CONJG(vme(1,ivbm,ivbm+1,i3))*vme(1,ivbm,ivbm+1,i3)),REAL(CONJG(vme(2,ivbm,ivbm+1,i3))*vme(2,ivbm,ivbm+1,i3)) 
     !PRINT*, ( REAL(CONJG(vme_LS(1,2*ivbm,2*ivbm+1,i3))*vme_LS(1,2*ivbm,2*ivbm+1,i3)) + &
     ! & REAL(CONJG(vme_LS(1,2*ivbm-1,2*ivbm+1,i3))*vme_LS(1,2*ivbm-1,2*ivbm+1,i3)) +  &
     ! & REAL(CONJG(vme_LS(1,2*ivbm-1,2*ivbm+2,i3))*vme_LS(1,2*ivbm-1,2*ivbm+2,i3)) + & 
     !& REAL(CONJG(vme_LS(1,2*ivbm  ,2*ivbm+2,i3))*vme_LS(1,2*ivbm  ,2*ivbm+2,i3)) ) / 2.0
     !PRINT*, ( REAL(CONJG(vme_LS(2,2*ivbm,2*ivbm+1,i3))*vme_LS(2,2*ivbm,2*ivbm+1,i3)) + &
     !& REAL(CONJG(vme_LS(2,2*ivbm-1,2*ivbm+1,i3))*vme_LS(2,2*ivbm-1,2*ivbm+1,i3)) +  &
     !& REAL(CONJG(vme_LS(2,2*ivbm-1,2*ivbm+2,i3))*vme_LS(2,2*ivbm-1,2*ivbm+2,i3)) + & 
     !& REAL(CONJG(vme_LS(2,2*ivbm  ,2*ivbm+2,i3))*vme_LS(2,2*ivbm  ,2*ivbm+2,i3)) ) / 2.0
     !STOP

     DO i = 1, nomega

        DO n = 1, 2*ivbm
        DO m = 2*ivbm+1, 2*nw
           DO j = 1, 3

              !!!!--- SPECIAL FOR PbO: wannierizarion was done for GW eigenvalues
           absorption_LS_local(j,i) = absorption_LS_local(j,i) + REAL(CONJG(vme_LS(j,n,m,i3))*vme_LS(j,n,m,i3))* &
                & EXP( -(Elda_out_local_LS(m,i3)-Elda_out_local_LS(n,i3) -omega(i))**2/eta**2 )/SQRT(pi*eta**2) 
           END DO
        END DO
        END DO

        DO n = 1, ivbm
        DO m = ivbm+1, nw
           DO j = 1, 3

              !!!!--- SPECIAL FOR PbO: wannierizarion was done for GW eigenvalues
           absorption_noLS_local(j,i) = absorption_noLS_local(j,i) + REAL(CONJG(vme(j,n,m,i3))*vme(j,n,m,i3))* &
                & EXP( -(Elda_out_local(m,i3)-Elda_out_local(n,i3) -omega(i))**2/eta**2 )/SQRT(pi*eta**2) 
           END DO
        END DO
        END DO

     END DO



  END DO
END DO iloop

#ifdef MPI                                                                                                                                                                  
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)                                                                                                                                    
   CALL MPI_ALLREDUCE( absorption_LS_local, absorption_LS, 3*nomega, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD, ierr)                                          
   CALL MPI_ALLREDUCE( absorption_noLS_local, absorption_noLS, 3*nomega, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD, ierr)                                          
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)                                                                                                                                    
#else                                                                                                                                                                       
   absorption_LS = absorption_LS_local
   absorption_noLS = absorption_noLS_local
#endif 


   !-- prefactors

   absorption_LS = absorption_LS / (Vcell*REAL(nkint(1)*nkint(2)*nkint(3),r8))
   absorption_noLS = absorption_noLS / (Vcell*REAL(nkint(1)*nkint(2)*nkint(3),r8))

   !-- calculate epsilon_2
   DO i = 1, nomega
      IF (omega(i) .GT. 1.0E-6) absorption_noLS(:,i) = absorption_noLS(:,i) * 8.0_r8*pi*pi/(omega(i)**2)

      !! BUF in next line: prefactor of 8 includes an excess factor of 2 for spin
      !IF (omega(i) .GT. 1.0E-6) absorption_LS(:,i) = absorption_LS(:,i) * 8.0_r8*pi*pi/(omega(i)**2)
      IF (omega(i) .GT. 1.0E-6) absorption_LS(:,i) = absorption_LS(:,i) * 4.0_r8*pi*pi/(omega(i)**2)
   END DO

   
   IF (ident .EQ. 0 ) THEN
      OPEN(unit=2,FILE='absorption_LS.dat')
      DO i = 1, nomega
      WRITE(2,*) omega(i)*ha2eV, absorption_LS(:,i)
      END DO
      CLOSE(2)
   END IF

   IF (ident .EQ. 0 ) THEN
      OPEN(unit=2,FILE='absorption_noLS.dat')
      DO i = 1, nomega
      WRITE(2,*) omega(i)*ha2eV, absorption_noLS(:,i)
      END DO
      CLOSE(2)
   END IF

#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_FINALIZE(ierr)
#endif

  STOP
END PROGRAM main
