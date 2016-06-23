! This program calculate the DOS of an finite square lattice
! using the chebyshev polynomial (CBP) expansion of the Green's
! function (GF) and calculating the GF directly using the 
! matrix inversion. Then two results will compare.
!====================================
! A. Ahmadi             April/6/2015
!=====================================
program sqrl_dos
 
  implicit none
  
  integer :: m_order, m1, m2, Mlat, Nlat, dim, ie, ne, im, nnz, in, s_num
  real*8 :: b_cof = 0.d0, e_max=0.d0, e_min=0.d0, e = 0.d0
  real*8 :: u_cheb = 0.0, den = 0.d0,  aux = 0.d0
  real*8,  parameter :: pi = 3.141592653, eta = 1.0d-8
  ! the Hamiltonian, identity and GF matrix 
  complex*16, allocatable, dimension(:,:) :: h_d, gf, tau_l, tau_r, &
       sigma_l, sigma_r, g_l, h_eff
  complex*16, allocatable, dimension(:) :: u_trace, h_values
  integer, allocatable, dimension(:) :: indx, jndx
  intrinsic dconjg
  Mlat = 50; Nlat = 1000; dim = mlat*nlat ; e_min = -4.0; e_max = 4.0; ne = 400
  m_order = 200; nnz = 2 * nlat * ( 2*mlat - 1 ); s_num = int(nnz/2)
  allocate( h_d(dim, dim), h_eff(dim,dim), gf(dim, dim) )
  allocate( tau_l(mlat, mlat*nlat), tau_r(mlat, dim), g_l(mlat,mlat) )
  allocate( sigma_l(dim,dim), sigma_r(dim,dim) )
  allocate( u_trace(0:m_order) )
  allocate( indx(nnz), jndx(nnz), h_values(nnz) )
  !=============================
  !    0-----0-----0-----0-----0-----2
  !  2 |   4 |   6 |   8 |  10 |  PBC   
  !  1 |   3 |   5 |   7 |   9 |  PBC
  !    0-----0-----0-----0-----0-----1
  !=============================




  !=========================
  !- constructing Hamiltonian
  !=========================
  in = 0 ! counting parameter
  vertical :do m1 = 0, nlat - 1 
     do m2 = 1, mlat - 1 
        in = in + 1
        indx(in) = m1*mlat + m2
        jndx(in) = m1*mlat + m2 +1
        h_values(in) = (-1.d0, 0.d0)
        h_d( m1*mlat + m2, m1*mlat + m2 +1 ) = (-1.d0, 0.d0) * &
             exp( (0.0, 1.0)*m1*b_cof  ) 
        !- considering the pierel coefficient for b-field
        indx(s_num + in) = m1*mlat + m2 +1
        jndx(s_num + in) = m1*mlat + m2
        h_values(s_num + in) = dconjg( h_values(in) )
        h_d( m1*mlat + m2 +1, m1*mlat + m2 ) = &
             dconjg(h_d( m1*mlat + m2, m1*mlat + m2 +1 ) )
     end do
  end do vertical

  horizintal: do m1 = 0, nlat - 2
     do m2 = 1, mlat
        in = in + 1
        indx(in) = m1*mlat + m2
        jndx(in) = (m1 + 1)*mlat  + m2
        h_values(in) = (-1.d0, 0.d0)
        h_d( m1*mlat + m2, (m1 + 1)*mlat  + m2 ) = (-1.d0, 0.d0) * &
             exp( (0.d0, 1.d0)*m2*b_cof  )
        indx(s_num + in) = (m1 + 1)*mlat  + m2
        jndx(s_num +  in) =  m1*mlat + m2
        h_values(s_num + in) = dconjg( h_values(in) )
        h_d( (m1 + 1)*mlat  + m2 , m1*mlat + m2) = &
             dconjg( h_d( m1*mlat + m2, (m1 + 1)*mlat  + m2 ) )
     end do
  end do horizintal
  !==============
  horizontal_PBC : do m2 = 1, mlat
     in = in + 1
     indx(in) = (nlat-1)*mlat + m2
     jndx(in) = m2
     h_values(in) = (-1.0, 0.d0)
     h_d( (nlat-1)*mlat + m2, m2 ) =(-1.d0, 0.d0) * &
          exp( (0.d0, 1.d0)*m2*b_cof  )
     indx(s_num + in) = m2
     jndx(s_num + in) = (nlat-1)*mlat + m2
     h_values(s_num + in) = dconjg( h_values(in) )
     h_d( m2,  (nlat-1)*mlat + m2 ) = dconjg ( h_d( (nlat-1)*mlat + m2, m2 ) )
  end do horizontal_PBC

!=====================================================================
  !====================================
  !- constructing left nd right contacts
  !====================================
  Do m1 = 1, Mlat
     tau_L(m1, m1) = (-1.d0, 0.d0)
     tau_R(m1, m1 + (Nlat - 1) * Mlat) = (-1.d0, 0.d0)
  End do
!=====================================================================  
  !- First : evaluating the DOS by calculating GF directly using
  !-  matrix inversion (LAPACK)
!=====================================================================
  !- loop over energy spectrum
!!$  energy_1 : do ie = 1, ne
!!$     e = ie*(e_max - e_min)/ne + e_min
!!$     call lead(e, g_l, mlat)
!!$     sigma_r = matmul(transpose(tau_l), matmul(g_l, tau_l))
!!$     sigma_l = matmul(transpose(tau_r), matmul(g_l, tau_r))
!!$     call inversion(nlat*mlat, e*(1.0,0.0)*iden - H_d  &
!!$          - sigma_l - sigma_r +(0.0, 1.0)*eta*iden , gf)
!!$     call inversion(nlat*mlat, e*(1.0,0.0)*iden - H_d + (0.0, 1.0)*eta*iden , gf)
!!$     call m_trace(nlat*mlat, gf, tr)
!!$     write(*,*) e, (-1.0/pi)*aimag(tr)
!!$3    format(f5.2,3x, f7.5)
!!$  end do energy_1

  
  !====================================================================
  !- Second : evaluating the DOS using the Chebyshev expansion of the GF
  !====================================================================
 
  ! call matrix_display(gf,dim)
!===========================================================================
!!$ Since the chebyshev terms of the Hamiltonina are completely independent of 
!!$ the energy, all terms can be calculated and stored in an array.
!===========================================================================
  
  h_d = (0.24d0, 0.d0) * h_d ! rescaling Hamiltonian for CBP expansion
  h_values = (0.24d0, 0.d0) * h_values
  !call cheb_generate(indx, jndx, h_values, u_trace, dim, m_order)
  call cheb_trace(indx, jndx, h_values, u_trace, dim, m_order)

  !- loop over energy
  energy_2 : do ie = 0, ne-1
     e = -1.d0 + ie*2.0/ne
     aux = 0.d0
     do im = 0, m_order-1
        u_cheb = ( 1.0/sin(acos(e)) ) * sin( (im+1) * acos(e) )
        aux = aux + real(u_trace(im)) * u_cheb * j_core(m_order,im)
     end do
     den = (2.0/pi) * sqrt(1 - e**2) * aux
     write(*,3) ie, e,  den/dim
3    format(i3, 1x, f5.2, 1x, f6.2)
  end do energy_2
  

  deallocate(H_d, gf, tau_l, tau_r, g_l, sigma_r, sigma_l,  u_trace )
  deallocate( indx, jndx, h_values )

  
contains
!======================================================================
subroutine cheb_trace(indx, jndx, h_values, u_trace, dim, m_order)
  use blas_sparse
  implicit none
  
  integer, intent(in) :: dim, m_order, indx(:), jndx(:)
  complex*16, intent(in) :: h_values(:)
  complex*16, dimension(dim) :: u_mat_pre, u_mat, u_mat_nex
  complex*16, intent(out), dimension(0:m_order) :: u_trace
  complex*16, dimension(dim) :: aux, rand_state
  real*8, dimension(dim) :: phase
  integer :: im, istat, a, ir, r_max 

  r_max = 30
  u_trace = (0.d0,0.d0)
  
  call zuscr_begin(dim,dim,a,istat) ! creating the handle of SM(Sparse Matrix)
  call uscr_insert_entries(a, h_values, jndx, indx, istat)! insert the entries at once
  call uscr_end(a,istat)! finish constructing 
  
  stoch_trace :  do ir = 0, r_max
     u_mat_pre = (0.d0, 0.d0)
     call random_number(phase)
     rand_state = exp( (0.d0,1.d0) * phase * 2 * pi )
     u_mat = rand_state
     CBP_generator : do im = 0, m_order-1
        aux = (0.d0,0.d0)
        call usmv(a, u_mat, aux, istat) ! sparse matrix-vector multiplication
        u_trace(im) = u_trace(im) + dot_product(rand_state, u_mat)
        u_mat_nex = 2 * aux - u_mat_pre
        u_mat_pre = u_mat
        u_mat = u_mat_nex
        
     end do CBP_generator
  end do stoch_trace
  call usds(a, istat)
  u_trace = u_trace / r_max
end subroutine cheb_trace
!=======================================
!- subroutine to expand the Hamiltoninan
!- in terms of the Chebyshev Polynomial
!- The trace of CBP will be stored in u_trace
!=======================================
subroutine cheb_generate(indx, jndx, h_values, u_trace,dim, m_order)
  use blas_sparse
  implicit none
  integer, intent(in) :: dim, m_order, indx(:), jndx(:)
  complex*16, intent(in) :: h_values(:)
  complex*16, dimension(dim,dim) :: u_mat_pre, u_mat, u_mat_nex
  complex*16, intent(out), dimension(0:m_order) :: u_trace
  complex*16, dimension(dim,dim) :: aux 
  integer :: im, istat, a

  u_trace = (0.d0,0.d0)
  !- create sparse BLAS handle

  call zuscr_begin(dim,dim,a,istat)
  
  !- insert entries

  call uscr_insert_entries(a, h_values, jndx, indx, istat)

  !- complete construction

  call uscr_end(a,istat)

  !-initializing for U_0(H)
  h_eff = h_d 
  u_mat_pre =  (0.d0,0.d0)
  forall (im = 1 : dim) u_mat(im,im) = (1.d0, 0.d0)
  
  do im = 0, m_order-1
     call m_trace( dim, u_mat, u_trace(im) )
     aux = (0.d0,0.d0)
     !call usmm(a,u_mat,aux, istat)
     !if (istat .ne. 0) write(*,*) "Something is wrong!" 
     aux = matmul(u_mat,h_eff)

     u_mat_nex = 2 * aux  - u_mat_pre
     u_mat_pre = u_mat
     u_mat = u_mat_nex
  end do

  call usds(a,istat)

end subroutine cheb_generate

!====================
!- the subroutine to calculate the green's 
!- function of the left and right leads
!==================== 
subroutine lead_scale(e, g_l, mlat)
  
  implicit none

  integer :: mlat ! the number of vertical sites in one slice  

  complex*16, dimension(mlat,mlat) :: h, u 
  complex*16, intent(out) :: g_l (:,:)
  
  double precision :: e
  complex*16, parameter :: eta = (0.d0,1.d0)*1.0d-8   
  
  integer :: i
  
  h = (0.d0,0.d0) !initializing

  u = (0.d0,0.d0) !initializing
  
  !====================
  !- constructing the hamiltonian and coupling matrices
  !====================
  do i = 1, mlat - 1

     h(i,i+1) = (-1.d0,0.d0)
     h(i+1,i) = dconjg(h(i,i+1))
 
  end do

  !===================
  do i = 1, mlat

     u(i,i) = (-1.d0,0.d0)
 
  end do
  !============================
  !- rescalling h and u for CBP
  !============================
  h = (0.25d0,0.d0)*h
  u = (0.25d0,0.d0)*u
  
  g_l = (0.d0,0.d0)
  
  !============================================
  !- loop to calculate the lead green's function
  !=============================================
 
  call glead_dec(e,h,u,g_l,mlat)
  

end subroutine lead_scale
!==============================================
  subroutine matrix_display(A, Mlat)
    
    implicit none
    integer, intent(in) :: Mlat ! dimension of the square A matrix
    integer :: m1, m2 ! auxilary counter
    complex*16, dimension(Mlat, Mlat), intent(in) :: A
    
    Do m1 = 1, Mlat 
       write(*,2) (abs(A(m1,m2)), m2 = 1, Mlat)
2      format(40(f4.2, 2X))
    end Do
  end subroutine matrix_display
  !=============================================
  !matrix inversion subroutine using lapack library
  subroutine inversion(idim,a,ainv)
    
    
    implicit none
    
  !===========
    ! variables
    !===========
    ! idim : (input) dimension of square complex matrix
    ! a    : (input) target matrix that we want to find the invese
    ! ainv : (output) inverse of matrix a*ainv = i
    integer, intent(in) :: idim
  complex*16, intent(in), dimension(idim,idim) :: a
  complex*16, intent(out), dimension(idim,idim) :: ainv
  
  
  !===============
  !local variables
  !===============
  
  integer  :: info
  integer, dimension(idim) :: ipiv
  
  complex*16, dimension(idim,idim) :: b 
  complex*16, dimension(idim) :: work
  
  b = a
  info = 0
  
  call zgetrf( idim, idim, b, idim, ipiv, info)
  call zgetri( idim, b, idim, ipiv, work, idim, info)
  if (info == 0) then
     ainv = b
  else
     write(*,*) "error: matrix inversion"
     stop
  end if
end subroutine inversion
!==============================================================
! Claculating the trace of matrix
subroutine m_trace(dim, a_matrix, tr)
  
  implicit none
  
  integer, intent(in) :: dim ! dimension of the matrix
  complex*16, intent(in), dimension(dim,dim) :: a_matrix
  complex*16, intent(out) :: tr 
  integer :: i !dummy
  tr = (0.d0, 0.d0)
  do i = 1, dim
     tr = tr + a_matrix(i,i)
  end do
end subroutine m_trace
!===============================================================
! Calculating the Jackson core to prevent the Gibbs oscillation in summation
pure function j_core(max,m)
  integer, intent(in) :: max, m
  real*8 :: j_core, aux1, aux2, aux3
  aux1 = (max - m + 1) * cos( pi * m / ( max+1 ) )
  aux2 = sin( pi * m /(max+1) )
  aux3 = ( 1.0 / tan( pi/(max+1) ) ) 
  j_core = ( 1.0/(max +1) ) * ( aux1 + aux2 * aux3)
end function j_core
!================================================================
!!$subroutine sp_mult
!!$use blas_sparse
!!$
!!$implicit none
!!$integer nmax, nnz
!!$parameter (nmax = 4, nnz = 6)
!!$integer i, n, a, istat
!!$integer, dimension(:), allocatable :: indx, jndx
!!$double precision, dimension(:), allocatable :: val, x, y
!!$
!!$allocate( val(nnz), x(nmax), y(nmax), indx(nnz), jndx(nnz) )
!!$
!!$x = (/1,2,3,4/)
!!$y = 0
!!$indx = (/1,2,2,3,4,4/)
!!$jndx = (/1,2,4,3,1,4/)
!!$val = (/1.1,2.2,2.4,3.3,4.1,4.4/)
!!$
!!$n = nmax
!!$
!!$!   ----------------------------------
!!$!-  step 1 : Create Sparse BLAS Handle
!!$!   ---------------------------------
!!$
!!$call duscr_begin(n,n,a,istat)
!!$
!!$!   --------------------------------
!!$!-   step 2 : Insert entries one-by-one
!!$!   --------------------------------
!!$
!!$do i = 1, nnz
!!$   call uscr_insert_entry(A, val(i), indx(i), jndx(i), istat)
!!$end do
!!$
!!$!   -----------------------------
!!$!-    step 3: complete construction of sparse matrix
!!$!   -----------------------------
!!$
!!$call uscr_end(a,istat)
!!$
!!$!    -----------------------------------------
!!$!-  step 4 : Compute Matrix vector product y = A*x
!!$!    -----------------------------------------
!!$
!!$call usmv(a,x,y,istat)
!!$
!!$write(*,3) (y(i), i = 1, 4)
!!$3 format(10(f5.2,2x))
!!$
!!$!   --------------------------------
!!$!-  step 5 : Release Matrix Handle
!!$!   -------------------------------
!!$
!!$call usds(a,istat)
!!$
!!$
!!$
!!$end program f95_sparse_ex
!!$
!================================================================
!==================================
SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  
  CALL SYSTEM_CLOCK(COUNT=clock)
  
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed


include "Glead_dec.f95"
include  "Lead.f95"
end program sqrl_dos

