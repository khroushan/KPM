!=========================================================================
! This program calculate the conductance of an finite square lattice
! using the chebyshev polynomial (CBP) expansion of the Green's
! function (GF).
!====================================
! A. Ahmadi             April/6/2015
! A. Ahmadi             April/27/2015 
! A. Ahmadi             May/4/2015     using Kubo-Bastin formulation
! A. Ahmadi             July/19/2015   fixed the Peierl's coefficient
! A. Ahmadi             July/31/2015   fixed the Peierl's coefficient (probably)
! A. Ahmadi             Aug/5/2015     magnetic field with PBC
! A. Ahmadi             Aug/6/2015     the velocity got fixed
!=====================================
module numeric_kind
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: sp = kind(1.0), &
                        dp = selected_real_kind(2 * precision(1.0_sp) ), &
                        qp = selected_real_kind(2 * precision(1.0_dp) )
  real(kind=dp), parameter :: pi = 3.14159265358979323
contains
  !================================================================
  SUBROUTINE init_random_seed(rn_seed)
    integer, intent(in) :: rn_seed
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    
    CALL SYSTEM_CLOCK(COUNT=clock)
    
    seed = clock + rn_seed * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    
    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed
  
end module numeric_kind
!=========================================================================
module hamiltonian
  use numeric_kind
contains

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

  mag_coef = (0.d0,1.d0) * (xy_j(1) - xy_i(1)) * ( xy_j(2) + xy_i(2) )

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
  character(len=30) :: fmt
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
        hxxh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(1) - xy_j(1) )
        hyyh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(2) - xy_j(2) )

        indx(s_num + in) = m1*mlat + m2 +1
        jndx(s_num + in) = m1*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
        hxxh_v(s_num + in) = conjg( hxxh_v(in) ) 
        hyyh_v(s_num + in) = conjg( hyyh_v(in) ) 
     end do
  end do vertical
  !==============
  vertical_PBC : do m1 = 1, nlat
     in = in + 1
     indx(in) = m1 * mlat
     jndx(in) = (m1-1) * mlat + 1
     xy_i(1) = x_cor(m1,mlat); xy_i(2) = y_cor(m1,mlat) ! coordinates of sites i
     xy_j(1) = x_cor(m1,1); xy_j(2) = y_cor(m1,1) ! coordinates of sites j
     h_values(in) = -(1.d0,0.d0) 
     hxxh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(1) - xy_j(1) )
     hyyh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(2) - xy_j(2) )
     
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
           hxxh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(1) - xy_j(1) ) !- construct [H,x], [H,y] 
           hyyh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(2) - xy_j(2) )
           
           indx(s_num + in) = (m1 + 1)*mlat  + ( m2 + 1 ) 
           jndx(s_num +  in) =  m1*mlat + m2
           h_values(s_num + in) = conjg( h_values(in) )
           hxxh_v(s_num + in) = conjg( hxxh_v(in) )
           hyyh_v(s_num + in) = conjg( hyyh_v(in) )
        else if ( mod(m2,4) == 0 ) then
           in = in + 1
           indx(in) = m1*mlat + m2
           jndx(in) = (m1 + 1)*mlat  + ( m2 - 1 )
           xy_i(1) = x_cor(m1+1,m2); xy_i(2) = y_cor(m1+1,m2) ! coordinates of sites i
           xy_j(1) = x_cor(m1+2,m2-1); xy_j(2) = y_cor(m1+2,m2-1) 
           h_values(in) = -exp( flux * mag_coef(xy_i, xy_j) )
           hxxh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(1) - xy_j(1) )
           hyyh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(2) - xy_j(2) )
           
           indx(s_num + in) = (m1 + 1)*mlat  + ( m2 - 1 ) 
           jndx(s_num +  in) =  m1*mlat + m2
           h_values(s_num + in) = conjg( h_values(in) )
           hxxh_v(s_num + in) = conjg( hxxh_v(in) )
           hyyh_v(s_num + in) = conjg( hyyh_v(in) )
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
        hxxh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(1) - xy_j(1) )
        hyyh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(2) - xy_j(2) )
        
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
        hxxh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(1) - xy_j(1) )
        hyyh_v(in) = (0.d0,1.d0) * h_values(in) * ( xy_i(2) - xy_j(2) )

        
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

  write(fmt,'(a,i3,a)') "(",nnz,"(i2,1x))"
  Open(unit = 15, file = "hamiltonian.dat", action = 'write')
  write(15,*) nnz
  write(15,fmt) indx
  write(15,*)
  write(15,fmt) jndx
  close(unit=15)
  
 
  deallocate(x_cor, y_cor)
end subroutine hamilton_graphene

  !========================================================
  !- constructing the Hamiltonina for an square lattice 
  !========================================================
subroutine  hamilton_sqr(mlat, nlat, flux, imp_amp, rn_seed, indx, jndx, &
       h_values, hxxh_v, hyyh_v)

  implicit none
  integer, intent(in) :: mlat, nlat, rn_seed
  integer, intent(out) :: indx(:), jndx(:)
  complex(kind=dp), intent(out) :: h_values(:), hxxh_v(:), hyyh_v(:)
  real(kind=dp), intent(in) :: flux, imp_amp 
  real(kind=dp) ::  rand_aux = 0.0 
  integer :: m1, m2, in, dim, s_num

  dim = mlat*nlat !dimension of the matrix
  s_num = 2 * dim

  !=============================
  !    0-----0-----0-----0-----0-----2
  !  2 |   4 |   6 |   8 |  10 |  PBC   
  !  1 |   3 |   5 |   7 |   9 |  PBC
  !    0-----0-----0-----0-----0-----1
  !=============================
  h_values = (0.d0, 0.d0); hxxh_v = (0.d0, 0.d0); hyyh_v = (0.d0, 0.d0)
  indx = 0; jndx = 0
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
        hyyh_v(in) = (0.d0,1.d0) * h_values(in)

        indx(s_num + in) = m1*mlat + m2 +1
        jndx(s_num + in) = m1*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
        hyyh_v(s_num + in) = conjg( hyyh_v(in) )
     end do
  end do vertical

!==============
  vertical_PBC : do m1 = 1, nlat
     in = in + 1
     indx(in) = m1 * mlat
     jndx(in) = (m1-1) * mlat + 1
     h_values(in) = (-1.0, 0.d0)
     hyyh_v(in) = (0.d0,1.d0) * h_values(in)

     indx(s_num + in) = (m1-1) * mlat + 1
     jndx(s_num + in) = m1 * mlat
     h_values(s_num + in) = conjg( h_values(in) )
     hyyh_v(s_num + in) = conjg( hyyh_v(in) )
  end do vertical_PBC
  !======================================================

  horizintal: do m1 = 0, nlat - 2
     do m2 = 1, mlat
        in = in + 1
        indx(in) = m1*mlat + m2
        jndx(in) = (m1 + 1)*mlat  + m2
        h_values(in) = (-1.d0, 0.d0) * exp( (0.d0,1.d0) * m2 * pi * flux )
        hxxh_v(in) = (0.d0,1.d0) * h_values(in)

        indx(s_num + in) = (m1 + 1)*mlat  + m2
        jndx(s_num + in) =  m1*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
        hxxh_v(s_num + in) = conjg( hxxh_v(in) )
        
     end do
  end do horizintal
  !==============

  horizontal_PBC : do m2 = 1, mlat
     in = in + 1
     indx(in) = (nlat-1)*mlat + m2
     jndx(in) = m2
     ! Not sure if the vector potential must be presented in the PBC
     h_values(in) = (-1.0, 0.d0) * exp( (0.d0,1.d0) * m2 * pi * flux )
     hxxh_v(in) = (0.d0,1.d0) * h_values(in)

     indx(s_num + in) = m2
     jndx(s_num + in) = (nlat-1)*mlat + m2
     h_values(s_num + in) = conjg( h_values(in) )
     hxxh_v(s_num + in) = conjg( hxxh_v(in) )
  end do horizontal_PBC
  
  call  init_random_seed(rn_seed)
  
  onsite : do m1 = 1, dim
     indx( 4*dim + m1 ) = m1
     jndx( 4*dim + m1 ) = m1
     call random_number(rand_aux)
     h_values( 4*dim + m1 ) = imp_amp * (rand_aux - 0.5)
  end do onsite

!!$  do in = 1, size(h_values)
!!$     hxxh_v(in) = h_values(in) * ( jndx(in) / mlat - indx(in) / mlat )
!!$     hyyh_v(in) = h_values(in) * ( mod(jndx(in),mlat) - mod(indx(in),mlat) )
!!$  end do

end subroutine hamilton_sqr

end module hamiltonian
!===========================================================================
module KPM_routines
  use numeric_kind
contains
  !================================================
!- Claculate the |n> using  CBP recursive formula
!- and store them in the array r_n
!=================================================
subroutine cheb_r(indx, jndx, h_values, r_n, m_order, r_0)
  implicit none
  
  integer, intent(in) :: m_order, indx(:), jndx(:)
  complex(kind=dp), intent(in) :: h_values(:), r_0(:)
  complex(kind=dp), intent(out) :: r_n(:,:)
  complex(kind=dp), dimension(size(r_0)) :: aux
  integer :: im
  

  r_n(1,:) = r_0
  call smvm(indx, jndx, h_values, r_n(1,:), aux)
  r_n(2,:) = 2 * aux
  CBP_generator : do im = 2, m_order - 1
     call smvm(indx, jndx, h_values, r_n(im,:), aux)
     r_n(im + 1,:) = 2 * aux - r_n(im-1,:)
  end do CBP_generator
end subroutine cheb_r
!=============================================================
!- Calculate the trace of gamma function which is defined sa
!-    gamma = Tr[hbar v_a * U_m(H) * hbar v_b * U_n(H) ] 
!=============================================================
subroutine gamma( indx, jndx, h_values, hxxh_v, hyyh_v, gamma_trace, &
     dim, m_order, r_max)
  implicit none
  
  integer, intent(in) :: r_max, dim, m_order, indx(:), jndx(:)
  complex(kind=dp), intent(in) :: h_values(:), hxxh_v(:), hyyh_v(:)  
  complex(kind=dp), dimension(m_order, dim) :: r_n
  complex(kind=dp), intent(out) :: gamma_trace(:,:)
  complex(kind=dp), dimension(dim) :: aux, rbar_n, rbar_np1, &
       rbar_nm1, r_0
  real(kind=dp), dimension(dim) :: phase
  integer :: im, ir, in 
  
  gamma_trace = (0.d0,0.d0)
  

  stoch_trace :  do ir = 1, r_max
     
     call random_number(phase)
     r_0 = exp( (0.d0,1.d0) * phase * 2 * pi )
     
     call smvm(indx, jndx, hxxh_v, r_0, rbar_n)
     
     !- generate and store |r >_n
     call cheb_r(indx, jndx, h_values, r_n, m_order, r_0)
     !-generate rbar and dot with |r>
     rbar_nm1 = (0.d0, 0.d0)
     CBP_recurssion : do im = 1, m_order
        
        do in = 1, m_order 
           call smvm(indx, jndx, hyyh_v, r_n(in,:), aux)
           gamma_trace(im,in) = gamma_trace(im,in) + dot_product(rbar_n, aux)   
        end do

        aux = (0.d0,0.d0)
        call smvm(indx, jndx, h_values, rbar_n, aux)
        rbar_np1 = 2 * aux - rbar_nm1
        rbar_nm1 = rbar_n
        rbar_n = rbar_np1
     end do CBP_recurssion
  end do stoch_trace
  gamma_trace = gamma_trace / r_max
end subroutine gamma



!=====================================
!- Sparse matrix-vector multiplication
!=====================================
subroutine smvm(indx, jndx, values, vector, vec_reslt)
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
!==============================================
subroutine matrix_display(A, Mlat)
    
  implicit none
  integer, intent(in) :: Mlat ! dimension of the square A matrix
  integer :: m1, m2 ! auxilary counter
  complex(kind=dp), dimension(Mlat, Mlat), intent(in) :: A
  
  Do m1 = 1, Mlat 
     write(*,2) (abs(A(m1,m2)), m2 = 1, Mlat)
2    format(40(f4.2, 2X))
  end Do
end subroutine matrix_display

!===============================================================
! Calculating the Jackson core to prevent the Gibbs oscillation 
!===============================================================
pure function j_core(max,m)
  integer, intent(in) :: max, m
  real(kind=dp) :: j_core, aux1, aux2, aux3
  aux1 = (max - m + 1) * cos( pi * m / ( max+1 ) )
  aux2 = sin( pi * m /(max+1) )
  aux3 = ( 1.0 / tan( pi/(max+1) ) ) 
  j_core = ( 1.0/(max +1) ) * ( aux1 + aux2 * aux3)
end function j_core

end module KPM_routines

!=========================================================================
! This program calculate the conductance of an finite square lattice
! using the chebyshev polynomial (CBP) expansion of the Green's
! function (GF).
!====================================
! A. Ahmadi             April/6/2015
! A. Ahmadi             April/27/2015 
! A. Ahmadi             MAy/4/2015     using Kubo-Bastin formulation
! A. Ahmadi             July/19/2015   fixed the Peierl's coefficient
! A. Ahmadi             July/31/2015   fixed the Peierl's coefficient (probably)
! A. Ahmadi             Aug/5/2015     magnetic field with PBC
!=====================================
program sqrl_conduct
 
  use numeric_kind
  use hamiltonian
  use kpm_routines
  implicit none
  
  
  integer :: m_order, mlat, nlat, dim, ie, ne, im, nnz, i_lattice = 0, &
       in, r_max, rn_seed = 37
  real(kind=dp) :: u_cheb_m = 0.0, u_cheb_n = 0.0, flux = 0.d0, imp_amp,  e = 0.d0
  ! the Hamiltonian, identity and GF matrix 
  complex(kind=dp) :: aux = (0.d0,0.d0), aux1 = (0.d0,0.d0), conduct = (0.d0, 0.d0)
  complex(kind=dp), allocatable, dimension(:) :: h_values, hxxh_v, hyyh_v
  complex(kind=dp), allocatable, dimension(:,:) :: gamma_trace
  integer, allocatable, dimension(:) :: indx, jndx
  character(len=20) :: file_in, file_out 
  
  if ( command_argument_count() == 1 ) then
     call get_command_argument(1, file_out)
     file_in = "CON_KPM_input.dat"
  elseif ( command_argument_count() == 2 ) then
     call get_command_argument(1, file_in)
     call get_command_argument(2, file_out)
  else
     file_in = 'CON_KPM_input.dat'
     file_out = 'result.dat'
  end if
  
  open(unit = 1, file = file_in , action = 'read')
  read(1,*) i_lattice        ! 0 = square, 1 = graphene
  read(1,*) mlat, nlat       ! number of vertical and horizontal sites    
  read(1,*) ne               ! energy resolution                         
  read(1,*) m_order          ! number of CBP to expand in terms of them  
  read(1,*) r_max            ! number of random vectors to evaluate trace
  read(1,*) flux             ! magnetic flux in unit of quantum flux
  read(1,*) imp_amp          ! amplitude of the anderson impurities
  read(1,*) rn_seed          ! external random seed 
  close(unit = 1)

  dim = mlat * nlat
  ! total number of nonzero elements in matrix (double PBC + onsite) 
  if (i_lattice == 0) then
     nnz =  5 * mlat * nlat
  elseif (i_lattice == 1) then
     if (mod(mlat,4) /= 0) then
        print *, "error mlat % 4 must be 0"
        stop
     else
        nnz = 4 * mlat * nlat
     end if
  else
     print *, " Error , 0 or 1....!"
     stop
  end if
  
  allocate( gamma_trace(m_order, m_order) )
  allocate( indx(nnz), jndx(nnz), h_values(nnz), hxxh_v(nnz), hyyh_v(nnz) )
  if (i_lattice == 0) then
     call hamilton_sqr(mlat, nlat, flux, imp_amp, rn_seed, indx, jndx, &
          h_values, hxxh_v, hyyh_v)
  elseif (i_lattice == 1) then
     call hamilton_graphene(mlat, nlat, flux, imp_amp, rn_seed, indx, jndx, &
          h_values, hxxh_v, hyyh_v )
  end if
  !====================================================================
  !- Evaluating the DOS using the Chebyshev expansion of the GF
  !====================================================================
  !rescalling the hamiltonian and velocities
  h_values = (0.25d0, 0.d0) * h_values
  hxxh_v = (0.25d0, 0.d0) * hxxh_v
  hyyh_v = (0.25d0, 0.d0) * hyyh_v

  !- calculate gamma_{mn} and store them in gamma_trace
  call gamma(indx, jndx, h_values, hxxh_v, hyyh_v, gamma_trace, dim, &
       m_order, r_max)
  
  deallocate( indx, jndx, h_values, hxxh_v, hyyh_v )
  conduct = (0.d0,0.d0)
  !- integration over  energy
  open(unit =2, file=file_out, action='write')
  energy_int : do ie = 0, ne
     e = pi * (ie + 0.5_dp) / ne
     !aux = 0.d0 because of integration over energy
     do im = 1, m_order
        u_cheb_m =  sin( im * e )
        do in = 1, m_order
           u_cheb_n =  sin( in * e )
           aux1 = in * u_cheb_m * exp( (0.d0,1.d0) * in * e ) + &
                  im * u_cheb_n * exp( (0.d0,-1.d0) * im * e )
           aux = aux +  aux1 * gamma_trace(im,in)  &
                * j_core(m_order,im - 1) * j_core(m_order,in - 1)
        end do
     end do
     conduct =  -8.0  * aux * (1.0/ne)
     write(2,3) ie, cos(e), real(conduct)/dim, aimag(conduct)/dim
3    format(i3, 1x, f6.3, 1x, f9.4, 1x, f9.4)
  end do energy_int
  close(unit=2)
  
  
  deallocate( gamma_trace )
  

end program sqrl_conduct

