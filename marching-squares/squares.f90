program marching_squares

    use glob
    implicit none
    integer :: i, j, k, k1, k2, k3, signcode
    integer, dimension(:), allocatable :: signs !vertex sign wrt. isovalue at each gridpoint

    call gridspecInput(xmin, xmax, dx, m, ymin, ymax, dy, n, sigma, xgrid, ygrid, fgrid) !user input to get init data
    call dataGen(xgrid, m, ygrid, n, fgrid) !generates function values at gridpoints
    call writeToFile(xmin, dx, m, ymin, dy, n) !writes info to 'grid_spec.txt' for plotting
    allocate(signs(m*n))
    
    do k = 1, m*n !calculate signs of each vertex
        if (fgrid(k) >= sigma) then
            signs(k) = 1;
        else 
            signs(k) = 0;
        end if 
    end do

    do i = 1, n !calculate interpolants of horizontal edges
        do j = 1, m-1
            call Index1D(k, m, i, j)
            call Index1D(k1, m-1, i, j)

            call LinearInterpolation(xgrid(j), fgrid(k), xgrid(j+1), fgrid(k+1), interpol_h(k1), sigma)

        end do
    end do 

    do i = 1, n-1 !calculate interpolants of vertical edges
        do j = 1, m
            call Index1D(k, m, i, j)
            call Index1D(k1, m, i+1, j)

            call LinearInterpolation(ygrid(i), fgrid(k), ygrid(i+1), fgrid(k1), interpol_v(k), sigma)
            
        end do
    end do 

    open (unit = 1, file = "edges.txt", status = "replace")

    do i = 1, n-1 ! loop to identify square type and construct edges
        do j = 1, m-1
            call Index1D(k1, m, i+1, j)
            call Index1D(k2, m, i+1, j+1)
            call Index1D(k3, m, i, j+1)
            call Index1D(k, m, i, j)

            signcode = 1 + signs(k1) + 2*signs(k2) + 4*signs(k3) + 8*signs(k); !converts to number from 1-16
            call constructEdges(signcode, i, j, xgrid(j), ygrid(i))
        end do
    end do

    close(1)

end program marching_squares

subroutine LinearInterpolation(x1, f1, x2, f2, x, f) !reverted linear interpolation
    implicit none
    real, intent(in) :: x1, f1, x2, f2, f
    real, intent(out) :: x

    x = x1 + (f-f1)*(x2-x1)/(f2-f1);

end subroutine LinearInterpolation

subroutine Index1D(k, n, i, j) !takes 2D index & #columns; returns 1D index by flattening row-wise
    implicit none
    integer, intent(in) :: i, j, n
    integer, intent(out) :: k

    k = n*(i-1) + j

end subroutine Index1D

subroutine Index2D(k, n, i, j) !reverse flatten, returns 2D index using 1D index 
    implicit none
    integer, intent(in) :: k, n
    integer, intent(out) :: i, j
    
    i = 1 + k/n
    j = modulo(k, n)

end subroutine Index2D

subroutine constructEdges(code, i, j, x, y)
    use glob
    implicit none
    integer, intent(in) :: code, i, j
    real, intent(in) :: x, y
    integer :: k1, k2, k3, k4
    character(len =10) :: fmt = "(*(1x,g0))"

    select case(code)
        case(2)
            call Index1D(k1, m, i, j)
            call Index1D(k2, m-1, i+1, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y+dy
        case(3)
            call Index1D(k1, m, i, j+1)
            call Index1D(k2, m-1, i+1, j)
            write(1,fmt) x+dx, interpol_v(k1), interpol_h(k2), y+dy
        case(4)
            call Index1D(k1, m, i, j)
            call Index1D(k2, m, i, j+1)
            write(1,fmt) x, interpol_v(k1), x+dx, interpol_v(k2)
        case(5)
            call Index1D(k1, m, i, j+1)
            call Index1D(k2, m-1, i, j)
            write(1,fmt) x+dx, interpol_v(k1), interpol_h(k2), y
        case(6)
            call Index1D(k1, m, i, j)
            call Index1D(k2, m-1, i, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y

            call Index1D(k3, m, i, j+1)
            call Index1D(k4, m-1, i+1, j)
            write(1,fmt) x+dx, interpol_v(k3), interpol_h(k4), y+dy
        case(7)
            call Index1D(k1, m-1, i, j)
            call Index1D(k2, m-1, i+1, j)
            write(1,fmt) interpol_h(k1), y, interpol_h(k2), y+dy
        case(8)
            call Index1D(k1, m, i, j)
            call Index1D(k2, m-1, i, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y
        case(9)
            call Index1D(k1, m, i, j)
            call Index1D(k2, m-1, i, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y
        case(10)
            call Index1D(k1, m-1, i, j)
            call Index1D(k2, m-1, i+1, j)
            write(1,fmt) interpol_h(k1), y, interpol_h(k2), y+dy
        case(11)
            call Index1D(k1, m, i, j)
            call Index1D(k2, m-1, i+1, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y+dy

            call Index1D(k3, m, i, j+1)
            call Index1D(k4, m-1, i, j)
            write(1,fmt) x+dx, interpol_v(k3), interpol_h(k4), y
        case(12)
            call Index1D(k1, m, i, j+1)
            call Index1D(k2, m-1, i, j)
            write(1,fmt) x+dx, interpol_v(k1), interpol_h(k2), y
        case(13)
            call Index1D(k1, m, i, j)
            call Index1D(k2, m, i, j+1)
            write(1,fmt) x, interpol_v(k1), x+dx, interpol_v(k2)
        case(14)
            call Index1D(k1, m, i, j+1)
            call Index1D(k2, m-1, i+1, j)
            write(1,fmt) x+dx, interpol_v(k1), interpol_h(k2), y+dy
        case(15)
            call Index1D(k1, m, i, j)
            call Index1D(k2, m-1, i+1, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y+dy
        
        case default

    end select


end subroutine constructEdges
