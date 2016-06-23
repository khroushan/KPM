module numeric_kind
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: sp = kind(1.0), &
                        dp = selected_real_kind(2 * precision(1.0_sp) ), &
                        qp = selected_real_kind(2 * precision(1.0_dp) )
  real(kind=dp), parameter :: pi = 4.0_dp * atan( 1.0_dp )
contains
  !============================
SUBROUTINE init_random_seed(rd_seed)
  !============================
  integer, intent(in) :: rd_seed
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  
  CALL SYSTEM_CLOCK(COUNT=clock)
  
  seed = clock + rd_seed * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed

end module numeric_kind
!=========================================================================
module hamiltonian
  use numeric_kind
contains
  !===========================================================
  !- constructing the hamiltoninan of an finite square lattice
  subroutine  sqr_hamilt(mlat, nlat, indx, jndx, h_values, onside_e, flux)
    !===========================================================
    implicit none
    integer, intent(in) :: mlat, nlat
    integer, intent(inout), dimension(:) :: indx, jndx
    real(kind=dp), intent(in) :: onside_e, flux
    complex(kind=dp), intent(inout), dimension(:) :: h_values
    real(kind=dp) :: aux
    integer :: nnz, s_num, m1, m2, in , dim
    nnz = 5 * nlat * mlat! total number of nonzero elements in matrix
    s_num = 2 * nlat * mlat
    dim = nlat * mlat 
    !=============================
    !        8               nlat*mlat
    !  4 0-----0-----0-----0-----0-----4
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

    !==============
  vertical_PBC : do m1 = 1, nlat
     in = in + 1
     indx(in) = m1 * mlat
     jndx(in) = (m1-1) * mlat + 1
     h_values(in) = (-1.0, 0.d0)

     indx(s_num + in) = (m1-1) * mlat + 1
     jndx(s_num + in) = m1 * mlat
     h_values(s_num + in) = conjg( h_values(in) )
  end do vertical_PBC
  !======================================================
    
    horizintal: do m1 = 0, nlat - 2
       do m2 = 1, mlat
          in = in + 1
          indx(in) = m1*mlat + m2
          jndx(in) = (m1 + 1)*mlat  + m2
          h_values(in) = (-1.d0, 0.d0) * exp( (0.d0,1.d0) * m2 * pi * flux )

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
    in = 0
    onsite : do m1 = 0, nlat - 1 ! centered potential LDOS
       do m2 = 1, mlat
          in = in + 1
          call random_number(aux)
          indx(4*dim + in) = m1 * mlat + m2
          jndx(4*dim + in) = m1 * mlat + m2
          if ( ( (m1 - nlat/2)**2 + (m2 - mlat/2)**2) < (mlat/10)**2 ) then
             h_values(4*dim + in) = onside_e * (-1.d-1, 0.d0) * aux
          else
             h_values(4*dim + in) = onside_e * (-1.d-1, 0.d0) * aux !(0.d0, 0.d0)
          end if
       end do
    end do onsite

  end subroutine sqr_hamilt
  !========================================================
!===============================================
!- generates the coordinates of carbon atoms 
!===============================================
subroutine coordinates(x_cor, y_cor)
  implicit none
  
  real(kind=dp), parameter :: a_lattice = 1.0d0, l_angle = pi/6.d0
  integer :: i, j
  real(kind=dp), intent(out) :: x_cor(:,:), y_cor(:,:)
  integer, dimension(2) :: loop 

  loop = shape(x_cor)
  
  do i = 1, loop(1)
     do j = 1, loop(2)

        if ( mod(j,4) == 1 ) then
           y_cor(i,j) = j/4 * 3 * a_lattice
           x_cor(i,j) = ( 2 * i + 1 ) * a_lattice * cos(l_angle)
        elseif ( mod(j,4) == 2 ) then
           y_cor(i,j) =  ( j/4 * 3 + sin(l_angle) ) * a_lattice 
           x_cor(i,j) = 2 * i * a_lattice * cos(l_angle);
        elseif ( mod(j,4) == 3 ) then
           y_cor(i,j) = ( j/4 * 3 + 1 + sin(l_angle) ) * a_lattice 
           x_cor(i,j) = 2 * i * a_lattice * cos(l_angle)
        elseif ( mod(j,4) == 0 ) then
           y_cor(i,j) =  ( (j/4 -1) * 3 + 2 ) * a_lattice 
           x_cor(i,j) = ( 2 * i + 1 ) * a_lattice * cos(l_angle)
        end if
    
     end do
  end do

  open(unit=64, file="coordinates_gr.xyz", action='write')
  write(64,*) loop(1) * loop(2)
  write(64,*)
  do i = 1, loop(1) 
     do j = 1, loop(2)

        write(64,*) "C", x_cor(i,j), y_cor(i,j), "0.0"
     end do
  end do


end subroutine coordinates

!======================================================  
!- evaluate the Pierel's magnetic coefficient
!======================================================
complex(kind=dp) function mag_coef(xy_i,xy_j)

  implicit none
  real(kind=8), intent(in) :: xy_i(:), xy_j(:)
  
  mag_coef = (0.d0,1.d0) * ( xy_j(2) - xy_i(2) ) * ( xy_j(1) + xy_i(1) )

end function mag_coef


!=====================================================
!- generate the Hamiltoninan matrix in sparse form wth
!- index stored in indx and jndx arrays
!=====================================================
subroutine hamilton_graphene(mlat, nlat, flux, imp_amp, rn_seed, indx, jndx, &
     h_values, hxxh_v, hyyh_v )

  implicit none
  integer, intent(in) :: mlat, nlat, rn_seed
  real(kind=dp), intent(in) :: flux, imp_amp
  integer, intent(out) :: indx(:), jndx(:)
  complex(kind=dp), intent(out) :: h_values(:), hxxh_v(:), hyyh_v(:)
  real(kind=dp), allocatable, dimension(:,:) :: x_cor, y_cor
  real(kind=dp), dimension(2) :: xy_i = 0.d0, xy_j = 0.d0
  real(kind=dp) :: rand_aux
  integer :: m1, m2, s_num, in, nnz , dim

  !=========================
  !- constructing Hamiltonian
  !=========================
  !        8               nlat*mlat
  !  4 0-----0-----0-----0-----0
  !    |\    |\    |\    |\    |\     
  !    |  \7 |  \  |  \  |  \  |  \ PBC 
  !  3 0    \0    \0    \0    \0    \ 3
  !    |     |     |     |     |
  !    |   6 |     |     |     |
  !  2 0    /0    /0    /0    /0    /2
  !    |  /  |  /  |  /  |  /  |  / 
  !    |/  5 |/    |/    |/    |/   PBC 
  !  1 0-----0-----0-----0-----0 
  !=============================
  allocate( x_cor(nlat,mlat), y_cor(nlat,mlat) )
  call coordinates(x_cor, y_cor)
  dim = mlat * nlat
  nnz = 4 * mlat * nlat 
  s_num = int(3 * mlat * nlat / 2)
  h_values = (0.d0, 0.d0); hxxh_v = (0.d0,0.d0); hyyh_v = (0.d0,0.d0)
  in = 0 ! counting parameter
  vertical :do m1 = 0, nlat - 1 
     do m2 = 1, mlat - 1 
        in = in + 1
        indx(in) = m1*mlat + m2
        jndx(in) = m1*mlat + m2 +1
        xy_i(1) = x_cor(m1+1,m2); xy_i(2) = y_cor(m1+1,m2) ! coordinates of sites i
        xy_j(1) = x_cor(m1+1,m2+1); xy_j(2) = y_cor(m1+1,m2+1) ! coordinates of sites j
        h_values(in) = -exp( flux * mag_coef(xy_i, xy_j) )
        hxxh_v(in) = h_values(in) * ( xy_i(1) - xy_j(1) )
        hyyh_v(in) = h_values(in) * ( xy_i(2) - xy_j(2) )

        indx(s_num + in) = m1*mlat + m2 +1
        jndx(s_num + in) = m1*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
        hxxh_v(s_num + in) = conjg(hxxh_v(in)) 
        hyyh_v(s_num + in) = conjg(hyyh_v(in)) 
     end do
  end do vertical
  !==============
  vertical_PBC : do m1 = 1, nlat
     in = in + 1
     indx(in) = m1 * mlat
     jndx(in) = (m1-1) * mlat + 1
     xy_i(1) = x_cor(m1,mlat); xy_i(2) = y_cor(m1,mlat) ! coordinates of sites i
     xy_j(1) = x_cor(m1,1); xy_j(2) = y_cor(m1,1) ! coordinates of sites j
     h_values(in) = -exp( flux * mag_coef(xy_i, xy_j) )!-(1.d0,0.d0) 
     hxxh_v(in) = h_values(in) * ( xy_i(1) - xy_j(1) )
     hyyh_v(in) = h_values(in) * ( xy_i(2) - xy_j(2) )
     
     indx(s_num + in) = (m1-1) * mlat + 1
     jndx(s_num + in) = m1 * mlat
     h_values(s_num + in) = conjg( h_values(in) )
     hxxh_v(s_num + in) = conjg( hxxh_v(in) )
     hyyh_v(s_num + in) = conjg( hyyh_v(in) )
  end do vertical_PBC
  !======================================================


  horizintal: do m1 = 0, nlat - 2
     do m2 = 1, mlat
        if ( mod(m2,4) == 1 ) then
           in = in + 1
           indx(in) = m1*mlat + m2
           jndx(in) = (m1 + 1)*mlat  + ( m2 + 1 )
           xy_i(1) = x_cor(m1+1,m2); xy_i(2) = y_cor(m1+1,m2) ! coordinates of sites i
           xy_j(1) = x_cor(m1+2,m2+1); xy_j(2) = y_cor(m1+2,m2+1) !coordinates of js
           h_values(in) = -exp( flux * mag_coef(xy_i, xy_j) )
           hxxh_v(in) = h_values(in) * ( xy_i(1) - xy_j(1) ) !- construct [H,x], [H,y] 
           hyyh_v(in) = h_values(in) * ( xy_i(2) - xy_j(2) )
           
           indx(s_num + in) = (m1 + 1)*mlat  + ( m2 + 1 ) 
           jndx(s_num +  in) =  m1*mlat + m2
           h_values(s_num + in) = conjg( h_values(in) )
           hxxh_v(s_num + in) = h_values(s_num + in) * ( xy_i(1) - xy_j(1) )
           hyyh_v(s_num + in) = h_values(s_num + in) * ( xy_i(2) - xy_j(2) )
        else if ( mod(m2,4) == 0 ) then
           in = in + 1
           indx(in) = m1*mlat + m2
           jndx(in) = (m1 + 1)*mlat  + ( m2 - 1 )
           xy_i(1) = x_cor(m1+1,m2); xy_i(2) = y_cor(m1+1,m2) ! coordinates of sites i
           xy_j(1) = x_cor(m1+2,m2-1); xy_j(2) = y_cor(m1+2,m2-1) 
           h_values(in) = -exp( flux * mag_coef(xy_i, xy_j) )
           hxxh_v(in) = h_values(in) * ( xy_i(1) - xy_j(1) )
           hyyh_v(in) = h_values(in) * ( xy_i(2) - xy_j(2) )
           
           indx(s_num + in) = (m1 + 1)*mlat  + ( m2 - 1 ) 
           jndx(s_num +  in) =  m1*mlat + m2
           h_values(s_num + in) = conjg( h_values(in) )
           hxxh_v(s_num + in) = h_values(s_num + in) * ( xy_i(1) - xy_j(1) )
           hyyh_v(s_num + in) = h_values(s_num + in) * ( xy_i(2) - xy_j(2) )
        end if
     end do
  end do horizintal
  !==============
  horizontal_PBC : do m2 = 1, mlat
     if ( mod(m2,4) == 1 ) then
        in = in + 1
        indx(in) = (nlat-1)*mlat + m2
        jndx(in) = m2 + 1        
        xy_i(1) = x_cor(nlat,m2); xy_i(2) = y_cor(nlat,m2) ! coordinates of sites i
        xy_j(1) = x_cor(1,m2+1); xy_j(2) = y_cor(1,m2+1) ! coordinates of sites j
        h_values(in) = -(1.d0,0.d0) 
        hxxh_v(in) = h_values(in) * ( xy_i(1) - xy_j(1) )
        hyyh_v(in) = h_values(in) * ( xy_i(2) - xy_j(2) )
        
        indx(s_num + in) = m2 + 1 
        jndx(s_num + in) = (nlat-1)*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
        hxxh_v(s_num + in) = conjg( hxxh_v(in) )
        hyyh_v(s_num + in) = conjg( hyyh_v(in) )
     else if ( mod(m2,4) == 0) then
        in = in + 1 
        indx(in) = (nlat-1)*mlat + m2
        jndx(in) = m2 - 1
        xy_i(1) = x_cor(nlat,m2); xy_i(2) = y_cor(nlat,m2) ! coordinates of sites i
        xy_j(1) = x_cor(1,m2-1); xy_j(2) = y_cor(1,m2-1) ! coordinates of sites j
        h_values(in) = -(1.d0,0.d0) 
        hxxh_v(in) = h_values(in) * ( xy_i(1) - xy_j(1) )
        hyyh_v(in) = h_values(in) * ( xy_i(2) - xy_j(2) )

        
        indx(s_num + in) = m2 - 1 
        jndx(s_num + in) = (nlat-1)*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
        hxxh_v(s_num + in) = conjg( hxxh_v(in) )
        hyyh_v(s_num + in) = conjg( hyyh_v(in) )
     end if
  end do horizontal_PBC
  call  init_random_seed(rn_seed)
  
  onsite : do m1 = 1, dim
     indx( 3*dim + m1 ) = m1
     jndx( 3*dim + m1 ) = m1
     call random_number(rand_aux)
     h_values( 3*dim + m1 ) = imp_amp * (rand_aux - 0.5)
  end do onsite
  
 
  deallocate(x_cor, y_cor)
end subroutine hamilton_graphene

  
end module hamiltonian
!=========================================================================

module matrix_cheb
  use numeric_kind

contains
  !=====================================
  !- Sparse matrix-vector multiplication
  !=====================================================
  subroutine smvm(indx, jndx, values, vector, vec_reslt)
  !=====================================================
    implicit none
    integer, intent(in) :: indx(:), jndx(:)
    complex(kind=dp), intent(in) :: values(:), vector(:)
    complex(kind=dp), intent(out) :: vec_reslt(:) 
    integer :: in
    vec_reslt = (0.d0, 0.d0)
    
    do in = 1, size(indx)
       vec_reslt( indx(in) ) = vec_reslt( indx(in) ) + values(in) * vector( jndx(in) )
    end do
  end subroutine smvm
  !===========================================
  !- sparse matrix-matrix calculation
  subroutine smmtrx(indx, jndx, a_values, b, c_res)
    integer, intent(in) :: indx(:), jndx(:) 
    complex(kind=dp), intent(in) :: a_values(:), b(:,:)
    complex(kind=dp), intent(out) :: c_res(:,:)
    integer :: i, in, dim_c(2)

    dim_c = shape(c_res)
    nnz = size(indx)
    c_res = 0.0
    
    do i = 1, dim_c(2)
       do in = 1, nnz
          c_res(indx(in),i) = c_res(indx(in),i) + a_values(in) * b(jndx(in),i)
       end do
    end do
  end subroutine smmtrx
  !===========================================================
  subroutine cheb_trace(indx, jndx, h_values, u_trace, dim, r_max)
  !===========================================================
    implicit none
  
    integer, intent(in) :: dim, r_max, indx(:), jndx(:)
    complex(kind=dp), intent(in) :: h_values(:)
    complex(kind=dp), dimension(dim) :: u_mat_pre, u_mat, u_mat_nex
    complex(kind=dp), intent(out) :: u_trace(:)
    complex(kind=dp), dimension(size(u_trace)) :: u_aux
    complex(kind=dp), dimension(dim) :: aux
    complex(kind=dp), dimension(r_max, dim) :: rand_state
    real(kind=dp), dimension(r_max,dim) :: phase
    integer :: im, ir, m_order
    
    m_order = size(u_trace)
    
    u_trace = (0.d0,0.d0); u_aux = (0.d0,0.d0)
    
    call random_number(phase)
    rand_state = exp( (0.d0,1.d0) * phase * 2 * pi )
    !$omp parallel do default(none), reduction(+ : u_trace), &
    !$omp& private(ir, im, u_mat, u_mat_pre, u_mat_nex, aux, u_aux), &
    !$omp& shared(r_max, m_order,rand_state, indx, jndx, h_values)
    stoch : do ir = 1, r_max   ! stochastic loop
       u_mat_pre = (0.d0, 0.d0)
       u_mat = rand_state(ir,:)
       CBP_generator : do im = 1, m_order
          aux = (0.d0,0.d0)
          call smvm(indx, jndx, h_values, u_mat, aux)
          u_aux(im) =  dot_product(rand_state(ir,:), u_mat)
          u_mat_nex = 2 * aux - u_mat_pre
          u_mat_pre = u_mat
          u_mat = u_mat_nex
       end do CBP_generator
       u_trace = u_trace + u_aux
    end do stoch
    !$omp end parallel do
    
    u_trace = u_trace / r_max
  end subroutine cheb_trace
  


  !===============================================================
  ! Calculating the Jackson core to prevent the Gibbs oscillation 
  !===============================================================
  pure function j_core(max,m) !
  !=============================
    integer, intent(in) :: max, m
    real(kind=dp) :: j_core, aux1, aux2, aux3
    aux1 = (max - m + 1) * cos( pi * m / ( max+1 ) )
    aux2 = sin( pi * m /(max+1) )
    aux3 = ( 1.0 / tan( pi/(max+1) ) ) 
    j_core = ( 1.0/(max +1) ) * ( aux1 + aux2 * aux3)
  end function j_core
  
end module matrix_cheb

!==============================================================
! This program calculate the DOS of an finite square lattice
! and graphene honeycomb lattice using the chebyshev polynomial 
! (CBP) expansion of the Green's function (GF)
!==============================================================
! A. Ahmadi             April/6/2015
! A. Ahmadi             June/12/2015         clean version using module
! A. Ahmadi             June/15/2015      parallel evaluation of cheb_trace, however it does not make it faster
! A. Ahmadi             June/24/2015      calculating the LDOS by eveluating the T(H) matrix recurssively
! A. Ahmai               Aug/5/2015       fixed the magnetic field in the graphene lattice
!=====================================
program ldos_kpm
  use numeric_kind
  use matrix_cheb
  use hamiltonian
!!  use omp_lib
  implicit none
  
  integer :: m_order, mlat, nlat, dim, ie, ne, im, in, nnz, r_max, cor
  real(kind=dp) :: u_cheb = 0.0, den = 0.d0,  aux = 0.d0, &
       e = 0.d0, e_ldos = 0.d0, onside_e = 0.d0, flux = 0.d0, imp_amp = 0.d0
  complex(kind=dp), allocatable, dimension(:) :: u_trace, h_values, hxxh_v, hyyh_v
  complex(kind=dp), allocatable, dimension(:,:) :: aux_mtx, u_mtx, u_mtx_p, u_mtx_n, rho_opt 
  real(kind=dp) , allocatable, dimension(:,:) :: x_cor, y_cor
  integer, allocatable, dimension(:) :: indx, jndx
  integer  :: lat_kind = 0, rn_seed = 37
  logical :: ldos_cal = .false. 
  character(len=15) :: file_name, fout_ldos
  intrinsic dconjg
  
  open(unit = 1, file = 'LDOS_input.dat', action = 'read')
  read(1,*) lat_kind         ! lattice kind ( 0 = square, 1 = graphene )
  read(1,*) mlat, nlat       ! number of horizental and verical sites    
  read(1,*) ne               ! energy resolution
  read(1,*) ldos_cal         ! ldos calculation (.true. = yes, .false. = no)
  read(1,*) e_ldos           ! energy of ldos
  read(1,*) onside_e         ! onside energy for few points
  read(1,*) flux             ! magnetic filed flux ( in unit of quantum flux )
  read(1,*) m_order          ! number of CBP to expand in terms of them  
  read(1,*) r_max            ! number of random vectors to evaluate trace
  read(1,*) file_name        ! name of file store dos using stochastic evaluation of trace
  read(1,*) fout_ldos        ! name of file to stor LDOS , x, y, ldos
  close(unit = 1)
  
  e_ldos = acos(e_ldos)
  dim = mlat * nlat ! dimension of matrix 
  lattice : if ( lat_kind == 0 ) then
     nnz = 5 * nlat * mlat ! total number of nonzero elements in matrix
     print *, "Square Lattice!"
  elseif ( lat_kind == 1 ) then
     nnz = 4 * mlat * nlat 
     print *, "Graphene lattice!"
  else
     print *, "error, the kind of lattice is unknown"
     stop  "valid values : (0) square (1) graphene"
  end if lattice
  
  allocate( u_trace(m_order) )
  allocate( indx(nnz), jndx(nnz), h_values(nnz), hxxh_v(nnz), hyyh_v(nnz) )
  
  !- constructing the Hamiltonian
  if ( lat_kind == 0 ) then
     call sqr_hamilt(mlat, nlat, indx, jndx, h_values, onside_e, flux)
  elseif ( lat_kind == 1 ) then
     !call graph_hamilt(mlat, nlat, indx, jndx, h_values)
     imp_amp = onside_e ! just to check
     call hamilton_graphene(mlat, nlat, flux, imp_amp, rn_seed, indx, jndx, &
          h_values, hxxh_v, hyyh_v )
  end if
  
  print *, "Matrix construction is done!"
  
  !=============================================================
  !- Evaluating the DOS using the Chebyshev expansion of the GF
  !=============================================================
 
  h_values = (0.25d0, 0.d0) * h_values !rescalling the hamiltonian
  
  call cheb_trace(indx, jndx, h_values, u_trace, dim, r_max)
  print *, " Gamma function calculation is done! "

  !- save input parameter in output file 
  open(unit=22, file = file_name, action='write')
  write(22,*)"#",  lat_kind,     " lattice kind ( 0 = square, 1 = graphene )"
  write(22,*)"#",  mlat, nlat,   " number of horizental and verical sites "
  write(22,*)'#',  ne,           " energy resolution"
  write(22,*)'#',  ldos_cal,     " ldos calculation (.true. = yes, .false. = no)"
  write(22,*)"#",  e_ldos,       " energy of ldos"
  write(22,*)"#",  onside_e,     " onside energy for few points"
  write(22,*)"#",  flux,         " magnetic filed flux ( in unit of quantum flux "
  write(22,*)"#",  m_order,      " number of CBP to expand in terms of them  "
  write(22,*)"#",  r_max,        " number of random vectors to evaluate trace"
  write(22,*)"#",  file_name,    " name of file store dos using stochastic evaluation of trace"
  write(22,*)"#",  fout_ldos,    " name of file to stor LDOS , x, y, ldos"
  write(22,*)"#################################################################"
  write(22,*)
  !====================================================================  
  !- loop over energy
  energy_2 : do ie = 0, ne - 1
     e = pi * ( ie + 0.5 ) / ne
     aux = 0.d0
     do im = 0, m_order-1
        u_cheb =  sin( (im+1) * e )
        aux = aux + real( u_trace(im+1) ) * u_cheb * j_core(m_order,im)
     end do
     den = (1.0/pi) *  aux
     write(22,3) ie, cos(e),  den/dim
3    format(i3, 1x, f5.2, 1x, f6.2)
  end do energy_2
  
  close(unit = 22)
  deallocate( u_trace )
 
  LDOS :if (ldos_cal .eqv. .true.) then ! LDOS calculation
     allocate( aux_mtx(dim,dim), u_mtx(dim,dim), u_mtx_p(dim,dim) )
     allocate( u_mtx_n(dim,dim), rho_opt(dim,dim) )
     allocate( x_cor(nlat,mlat), y_cor(nlat,mlat) )
     call coordinates(x_cor, y_cor)
     !- initializing the Chebychev matrices
     u_mtx_p = (0.d0,0.d0); u_mtx = (0.d0,0.d0)
     forall(im = 1:dim)
        u_mtx(im,im) = (1.d0,0.d0)
     end forall
     
     rho_opt = j_core(m_order,0) * sin( e_ldos ) * u_mtx
     do im = 1, m_order
        call smmtrx(indx, jndx, h_values, u_mtx, aux_mtx)
        u_mtx_n = 2 * aux_mtx - u_mtx_p
        u_mtx_p = u_mtx
        u_mtx = u_mtx_n
        u_cheb = sin( (im+1) * e_ldos )
        rho_opt = rho_opt + j_core(m_order,im) * u_cheb * u_mtx 
     end do
     rho_opt = 2/pi * rho_opt
     open(unit=41, file = fout_ldos, action='write')
     do im = 1, nlat
        do in = 1, mlat
           cor = (im - 1) * mlat + in
           write(41,*) x_cor(im,in), y_cor(im,in) , real( rho_opt(cor,cor) )
        end do
     end do
     close(unit=41)
     deallocate( aux_mtx, u_mtx, u_mtx_p, u_mtx_n, rho_opt )
  end if LDOS
  deallocate( indx, jndx, h_values )
  
end program ldos_kpm
