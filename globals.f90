module glob
    implicit none
    integer :: m, n !number of rows and columns of 2D grid
    real :: xmin, xmax, dx, ymin, ymax, dy, sigma !grid specs
    real, dimension(:), allocatable :: xgrid, ygrid, fgrid !grids
    real, dimension(:), allocatable :: interpol_v, interpol_h !grid of edge interpolations

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

    subroutine gridspecInput(xmin, xmax, dx, m, ymin, ymax, dy, n, sigma, xgrid, ygrid, fgrid) !input and assignment of global variables 
        real, intent(out) :: xmin, xmax, dx, ymin, ymax, dy, sigma
        real, intent(out), allocatable :: xgrid(:), ygrid(:), fgrid(:)
        integer, intent(out) :: m, n

        write(*,*) "xmin, xmax, m"; 
        read(*,*) xmin; read(*,*) xmax; read(*,*) m
        
        write(*,*) "ymin, ymax, n"; 
        read(*,*) ymin; read(*,*) ymax; read(*,*) n
        
        write(*,*) "isovalue:";
        read(*,*) sigma;

        allocate(xgrid(m))
        allocate(ygrid(n))
        allocate(fgrid(m*n))
        allocate(interpol_h((m-1)*n))
        allocate(interpol_v((n-1)*m))

        call Grid1D(xmin, xmax, m, dx, xgrid)
        call Grid1D(ymin, ymax, n, dy, ygrid)
      
    end subroutine gridspecInput

    subroutine dataGen(xgrid, m, ygrid, n, fgrid) !generate data for function value at gridpoints
        real, intent(in) :: xgrid(m), ygrid(n)
        real, intent(out) :: fgrid(m*n)
        integer, intent(in) :: m, n
        integer :: i, j

        do i = 1, n
            do j = 1, m
                fgrid(m*(i-1) + j) = xgrid(j)**4 - xgrid(j)**2 + ygrid(i)**2
            end do 
        end do
        
    end subroutine dataGen

    subroutine writeToFile(xmin, dx, m, ymin, dy, n) !test function for grid, circular contours
        real, intent(in) :: xmin, dx, ymin, dy
        integer, intent(in) :: m, n
        open (unit = 3, file = "grid_specs.txt", status = "replace")
        write(3,"(*(1x,g0))") xmin, dx, m
        write(3,"(*(1x,g0))") ymin, dy, n
        close(3)
    end subroutine writeToFile
end module glob