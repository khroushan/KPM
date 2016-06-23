! This program calculate the DOS of an finite square lattice
! using the chebyshev polynomial (CBP) expansion of the Green's
! function (GF) and calculating the GF directly using the 
! matrix inversion. Then two results will compare.
!====================================
! A. Ahmadi             April/6/2015
!=====================================
program sqrl_dos
 
  implicit none
  
  integer :: m_order, mlat, nlat, dim, ie, ne, im, nnz, r_max
  real*8 :: u_cheb = 0.0, den = 0.d0,  aux = 0.d0, e = 0.d0
  real*8,  parameter :: pi = 4 * atan(1.d0)
  complex*16, allocatable, dimension(:) :: u_trace, h_values
  integer, allocatable, dimension(:) :: indx, jndx
  intrinsic dconjg
  
  open(unit = 1, file = 'input.dat', action = 'read')

  read(1,*) mlat, nlat       ! number of horizental and verical sites    
  read(1,*) ne               ! energy resolution                         
  read(1,*) m_order          ! number of CBP to expand in terms of them  
  read(1,*) r_max            ! number of random vectors to evaluate trace

  close(unit = 1)

  dim = mlat*nlat !dimension of the matrix
  nnz = 2 * nlat * ( 2*mlat - 1 )! total number of nonzero elements in matrix
  
  allocate( u_trace(m_order) )
  allocate( indx(nnz), jndx(nnz), h_values(nnz) )

  !- constructing the Hamiltonian
  call sqr_hamilt(mlat, nlat, indx, jndx, h_values)
  
  !====================================================================
  !- Evaluating the DOS using the Chebyshev expansion of the GF
  !====================================================================
 
  h_values = (0.25d0, 0.d0) * h_values !rescalling the hamiltonian

  call cheb_trace(indx, jndx, h_values, u_trace, r_max)
  
  deallocate( indx, jndx, h_values )
  !- loop over energy
  
  energy_2 : do ie = 0, ne - 1
     e = pi * ( ie + 0.5 ) / ne
     aux = 0.d0
     do im = 0, m_order-1
        u_cheb =  sin( (im+1) * e )
        aux = aux + real( u_trace(im+1) ) * u_cheb * j_core(m_order,im)
     end do
     den = (2.0/pi) *  aux
     write(*,3) ie, cos(e),  den/dim
3    format(i3, 1x, f5.2, 1x, f6.2)
  end do energy_2
  

  deallocate( u_trace )

  
contains
!===========================================================
!- constructing the hamiltoninan of an finite square lattice
subroutine  sqr_hamilt(mlat, nlat, indx, jndx, h_values)
!===========================================================
  implicit none
  integer, intent(in) :: mlat, nlat
  integer, intent(inout), dimension(:) :: indx, jndx
  complex*16 , intent(inout), dimension(:) :: h_values
  integer :: nnz, s_num, m1, m2, in 
  nnz = 2 * nlat * ( 2*mlat - 1 )! total number of nonzero elements in matrix
  s_num = int(nnz/2)

  !=============================
  !        8               nlat*mlat
  !  4 0-----0-----0-----0-----0-----3
  !    |     |     |     |     |     
  !    |   7 |     |     |     |  PBC 
  !  3 0-----0-----0-----0-----0-----3
  !    |     |     |     |     |
  !    |   6 |     |     |     |
  !  2 0-----0-----0-----0-----0-----2
  !    |     |     |     |     |   
  !    |   5 |     |     |     |  PBC 
  !  1 0-----0-----0-----0-----0-----1
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
        
        indx(s_num + in) = m1*mlat + m2 +1
        jndx(s_num + in) = m1*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
     end do
  end do vertical
  
  horizintal: do m1 = 0, nlat - 2
     do m2 = 1, mlat
        in = in + 1
        indx(in) = m1*mlat + m2
        jndx(in) = (m1 + 1)*mlat  + m2
        h_values(in) = (-1.d0, 0.d0)
        indx(s_num + in) = (m1 + 1)*mlat  + m2
        jndx(s_num +  in) =  m1*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
     end do
  end do horizintal
  !==============
  horizontal_PBC : do m2 = 1, mlat
     in = in + 1
     indx(in) = (nlat-1)*mlat + m2
     jndx(in) = m2
     h_values(in) = (-1.0, 0.d0)
     
     indx(s_num + in) = m2
     jndx(s_num + in) = (nlat-1)*mlat + m2
     h_values(s_num + in) = conjg( h_values(in) )
  end do horizontal_PBC

end subroutine sqr_hamilt
!========================================================
subroutine cheb_trace(indx, jndx, h_values, u_trace, r_max)
  implicit none
  
  integer, intent(in) :: r_max, indx(:), jndx(:)
  complex*16, intent(in) :: h_values(:)
  complex*16, dimension(dim) :: u_mat_pre, u_mat, u_mat_nex
  complex*16, intent(out) :: u_trace(:)
  complex*16, dimension(dim) :: aux, rand_state
  real*8, dimension(dim) :: phase
  integer :: im, ir, m_order
  
  m_order = size(u_trace)

  u_trace = (0.d0,0.d0)
  
  stoch_trace :  do ir = 1, r_max
     u_mat_pre = (0.d0, 0.d0)
     call random_number(phase)
     rand_state = exp( (0.d0,1.d0) * phase * 2 * pi )
     u_mat = rand_state
     CBP_generator : do im = 1, m_order
        aux = (0.d0,0.d0)
        call smvm(indx, jndx, h_values, u_mat, aux)
        u_trace(im) = u_trace(im) + dot_product(rand_state, u_mat)
        u_mat_nex = 2 * aux - u_mat_pre
        u_mat_pre = u_mat
        u_mat = u_mat_nex
     end do CBP_generator
  end do stoch_trace
  u_trace = u_trace / r_max
end subroutine cheb_trace
!=====================================
!- Sparse matrix-vector multiplication
!=====================================
subroutine smvm(indx, jndx, values, vector, vec_reslt)
  implicit none
  integer, intent(in) :: indx(:), jndx(:)
  complex*16, intent(in) :: values(:), vector(:)
  complex*16, intent(out) :: vec_reslt(:) 
  integer :: in
  vec_reslt = (0.d0, 0.d0)

  do in = 1, size(indx)
     vec_reslt( indx(in) ) = vec_reslt( indx(in) ) + values(in) * vector( jndx(in) )
  end do
end subroutine smvm

!==============================================================
! Claculating the trace of matrix
!=============================================
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
! Calculating the Jackson core to prevent the Gibbs oscillation 
!===============================================================
pure function j_core(max,m)
  integer, intent(in) :: max, m
  real*8 :: j_core, aux1, aux2, aux3
  aux1 = (max - m + 1) * cos( pi * m / ( max+1 ) )
  aux2 = sin( pi * m /(max+1) )
  aux3 = ( 1.0 / tan( pi/(max+1) ) ) 
  j_core = ( 1.0/(max +1) ) * ( aux1 + aux2 * aux3)
end function j_core
!================================================================
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

end program sqrl_dos

