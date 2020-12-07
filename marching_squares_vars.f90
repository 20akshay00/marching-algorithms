module marching_variables
    implicit none
    integer :: m, n 
    real :: xmin, xmax, dx, ymin, ymax, dy
    real, dimension(:), allocatable :: xgrid, ygrid, fgrid

contains

    subroutine Grid1D(xmin, xmax, n, step, grid) !generates 1d grid
        implicit none
    
        real, intent(in) :: xmin, xmax
        integer, intent(in) :: n
        real, intent(out) :: grid(n)
        real, intent(out) :: step        
        integer :: i
    
        step = (xmax - xmin)/(n-1) !step size of grid
    
        do i = 1, n
        grid(i) = xmin + (i-1)*step;
        end do

    end subroutine Grid1D

    subroutine gridspecInput(xmin, xmax, dx, m, ymin, ymax, dy, n, xgrid, ygrid, fgrid) !input and assignment of global variables 
        real, intent(out) :: xmin, xmax, dx, ymin, ymax, dy
        real, intent(out), allocatable :: xgrid(:), ygrid(:), fgrid(:)
        integer, intent(out) :: m, n

        write(*,*) "xmin, xmax, m"; 
        read(*,*) xmin; read(*,*) xmax; read(*,*) m
        
        write(*,*) "ymin, ymax, n"; 
        read(*,*) ymin; read(*,*) ymax; read(*,*) n
        
        allocate(xgrid(m))
        allocate(ygrid(n))
        allocate(fgrid(m*n))

        call Grid1D(xmin, xmax, m, dx, xgrid)
        call Grid1D(ymin, ymax, n, dy, ygrid)
      
    end subroutine gridspecInput

    subroutine dataGen(xgrid, m, ygrid, n, fgrid) !generate data for function value at gridpoints
        real, intent(in) :: xgrid(m), ygrid(n)
        real, intent(out) :: fgrid(m*n)
        integer, intent(in) :: m, n
        integer :: i, j

        do i = 1,m
            do j = 1,n
                fgrid(m*(i-1) + j) = xgrid(i)**2 + ygrid(j)**2 
            end do 
        end do
        
    end subroutine dataGen

end module marching_variables