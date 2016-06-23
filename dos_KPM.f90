module numeric_kind
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: sp = kind(1.0), &
                        dp = selected_real_kind(2 * precision(1.0_sp) ), &
                        qp = selected_real_kind(2 * precision(1.0_dp) )
  real(kind=dp), parameter :: pi = 4.0_dp * atan( 1.0_dp )
end module numeric_kind
!=========================================================================
module hamiltonian
  use numeric_kind
contains
  !===========================================================
  !- constructing the hamiltoninan of an finite square lattice
  subroutine  sqr_hamilt(mlat, nlat, indx, jndx, h_values)
    !===========================================================
    implicit none
    integer, intent(in) :: mlat, nlat
    integer, intent(inout), dimension(:) :: indx, jndx
    complex(kind=dp), intent(inout), dimension(:) :: h_values
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
!=====================================================
!- generate the Hamiltoninan matrix in sparse form wth
!- index stored in indx and jndx arrays
!=====================================================
subroutine graph_hamilt(mlat, nlat, indx, jndx, h_values)
 
  implicit none
  integer, intent(in) :: mlat, nlat
  integer, intent(out) :: indx(:), jndx(:)
  complex(kind=dp), intent(out) :: h_values(:)
  integer :: m1, m2, in, nnz, s_num

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

  nnz = 3 * mlat * nlat 
  s_num = int(nnz/2)

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
        if (mod(m2,4) == 1) then
           in = in + 1
           indx(in) = m1*mlat + m2
           jndx(in) = (m1 + 1)*mlat  + ( m2 + 1 )
           h_values(in) = (-1.d0, 0.d0)
           
           indx(s_num + in) = (m1 + 1)*mlat  + ( m2 + 1 ) 
           jndx(s_num +  in) =  m1*mlat + m2
           h_values(s_num + in) = conjg( h_values(in) )
        else if (mod(m2,4) == 0 ) then
           in = in + 1
           indx(in) = m1*mlat + m2
           jndx(in) = (m1 + 1)*mlat  + ( m2 - 1 )
           h_values(in) = (-1.d0, 0.d0)
           
           indx(s_num + in) = (m1 + 1)*mlat  + ( m2 - 1 ) 
           jndx(s_num +  in) =  m1*mlat + m2
           h_values(s_num + in) = conjg( h_values(in) )
        end if
     end do
  end do horizintal
  !==============
  horizontal_PBC : do m2 = 1, mlat
     if ( mod(m2,4) == 1 ) then
        in = in + 1
        indx(in) = (nlat-1)*mlat + m2
        jndx(in) = m2 + 1
        h_values(in) = (-1.0, 0.d0) 
        
        indx(s_num + in) = m2 + 1 
        jndx(s_num + in) = (nlat-1)*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
     else if (mod(m2,4) == 0) then
        in = in + 1 
        indx(in) = (nlat-1)*mlat + m2
        jndx(in) = m2 - 1
        h_values(in) = (-1.0, 0.d0) 
        
        indx(s_num + in) = m2 - 1 
        jndx(s_num + in) = (nlat-1)*mlat + m2
        h_values(s_num + in) = conjg( h_values(in) )
     end if
  end do horizontal_PBC
  

end subroutine graph_hamilt  

  
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
    real(kind=dp), intent(in) :: a_values(:), b(:,:)
    real(kind=dp), intent(out) :: c_res(:,:)
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
  !============================
  SUBROUTINE init_random_seed()
  !============================
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    
    CALL SYSTEM_CLOCK(COUNT=clock)
    
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    
    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed

  
end module matrix_cheb

!==============================================================
! This program calculate the DOS of an finite square lattice
! and graphene honeycomb lattice using the chebyshev polynomial 
! (CBP) expansion of the Green's function (GF)
!==============================================================
! A. Ahmadi             April/6/2015
! A. Ahmadi             June/12/2015         clean version using module
! A. Ahmadi             June/15/2015      parallel evaluation of cheb_trace, however it does not make it faster
!=====================================
program sqrl_dos
  use numeric_kind
  use matrix_cheb
  use hamiltonian
!!  use omp_lib
  implicit none
  
  integer :: m_order, mlat, nlat, dim, ie, ne, im, nnz, r_max
  real(kind=dp) :: u_cheb = 0.0, den = 0.d0,  aux = 0.d0, e = 0.d0
  complex(kind=dp), allocatable, dimension(:) :: u_trace, h_values
  integer, allocatable, dimension(:) :: indx, jndx
  integer  :: lat_kind = 0
  character(len=15) :: file_name, ldos_file
  intrinsic dconjg
  
  open(unit = 1, file = 'KPM_input.dat', action = 'read')
  read(1,*) lat_kind         ! lattice kind ( 0 = square, 1 = graphene )
  read(1,*) mlat, nlat       ! number of horizental and verical sites    
  read(1,*) ne               ! energy resolution                         
  read(1,*) m_order          ! number of CBP to expand in terms of them  
  read(1,*) r_max            ! number of random vectors to evaluate trace
  read(1,*) file_name
  close(unit = 1)
  
  dim = mlat * nlat ! dimension of matrix 
  lattice : if ( lat_kind == 0 ) then
     nnz = 2 * nlat * ( 2*mlat - 1 )! total number of nonzero elements in matrix
     print *, "Square Lattice!"
  elseif ( lat_kind == 1 ) then
     nnz = 3 * mlat * nlat 
     print *, "Graphene lattice!"
  else
     print *, "error, the kind of lattice is unknown"
     stop  "valid values : (0) square (1) graphene"
  end if lattice
  
  allocate( u_trace(m_order) )
  allocate( indx(nnz), jndx(nnz), h_values(nnz) )
  
  !- constructing the Hamiltonian
  if ( lat_kind == 0 ) then
     call sqr_hamilt(mlat, nlat, indx, jndx, h_values)
  elseif ( lat_kind == 1 ) then
     call graph_hamilt(mlat, nlat, indx, jndx, h_values)
  end if
  
  
  !====================================================================
  !- Evaluating the DOS using the Chebyshev expansion of the GF
  !====================================================================
 
  h_values = (0.25d0, 0.d0) * h_values !rescalling the hamiltonian
  
  call cheb_trace(indx, jndx, h_values, u_trace, dim, r_max)
  
  deallocate( indx, jndx, h_values )
  !- save initial parameter on top of output file
  open(unit=22, file = file_name, action='write')
  write(22,*)"#", lat_kind,   "lattice kind ( 0 = square, 1 = graphene )"
  write(22,*)"#", mlat, nlat, "number of horizental and verical sites "
  write(22,*)"#", ne,         "        energy resolution  "
  write(22,*)"#", m_order,    "number of CBP to expand in terms of them  "
  write(22,*)"#", r_max,      "number of random vectors to evaluate trace"
  write(22,*) "###########################################################"
  write(22,*)
  !======================================================================
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
  
end program sqrl_dos
