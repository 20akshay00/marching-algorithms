program marching_squares

    use glob
    implicit none
    integer :: i, i_, j, k, k_, k1, k2, k3, signcode
    integer, dimension(:), allocatable :: signs !

    call gridspecInput(xmin, xmax, dx, m, ymin, ymax, dy, n, sigma, xgrid, ygrid, fgrid) !user input to get init data
    call dataGen(xgrid, m, ygrid, n, fgrid) !generates function values at gridpoints
    call writeToFile(xmin, dx, m, ymin, dy, n)
    allocate(signs(m*n))
    
    do k = 1, m*n !calculate signs of each vertex
        if (fgrid(k) >= sigma) then
            signs(k) = 1;
        else 
            signs(k) = 0;
        end if 
    end do

    do i = 1, m !calculate interpolants of horizontal edges
        do j = 1, n-1
            call Index1D(k, n, i, j)
            call Index1D(k_, n-1, i, j)

            call LinearInterpolation(xgrid(j), fgrid(k), xgrid(j+1), fgrid(k+1), interpol_h(k_), sigma)

        end do
    end do 

    do i = 1, m-1 !calculate interpolants of vertical edges
        i_ = i + 1
        do j = 1, n
            call Index1D(k, n, i, j)
            call Index1D(k_, n, i_, j)

            call LinearInterpolation(ygrid(i), fgrid(k), ygrid(i+1), fgrid(k_), interpol_v(k), sigma)
            
        end do
    end do 

    open (unit = 2, file = "grid.txt", status = "replace")
    write(2,"(*(1x,g0))") xgrid
    write(2,"(*(1x,g0))") fgrid
    write(2,"(*(1x,g0))") interpol_h
    write(2,"(*(1x,g0))") interpol_v
    close(2)

    open (unit = 1, file = "edges.txt", status = "replace")

    do i = 1, m-1 ! loop to identify square type and construct edges
        do j = 1, n-1
            call Index1D(k1, n, i+1, j)
            call Index1D(k2, n, i+1, j+1)
            call Index1D(k3, n, i, j+1)
            call Index1D(k, n, i, j)

            signcode = 1 + signs(k1) + 2*signs(k2) + 4*signs(k3) + 8*signs(k);
            call constructEdges(signcode, i, j, xgrid(k), ygrid(k))
        end do
    end do

    close(1)

end program marching_squares

subroutine LinearInterpolation(x1, f1, x2, f2, x, f) !reverted linear interpolation
    implicit none
    real, intent(in) :: x1, f1, x2, f2, f
    real, intent(out) :: x

    if((f1<f .and. f<f2) .or. (f2<f .and. f<f1)) then
        x = x1 + (f-f1)*(x2-x1)/(f2-f1);
    else
        x = 9000; !error code
    end if 

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
    integer, intent(in) :: code, i, j
    real, intent(in) :: x, y
    integer :: i1, j1, k1, k2, k3, k4, n1
    character(len =10) :: fmt = "(*(1x,g0))"

    i1 = i + 1
    j1 = j + 1
    n1 = n - 1

    select case(code)
        case(2)
            call Index1D(k1, n, i, j)
            call Index1D(k2, n1, i1, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y+dy
        case(3)
            call Index1D(k1, n, i, j1)
            call Index1D(k2, n1, i1, j)
            write(1,fmt) x+dx, interpol_v(k1), interpol_h(k2), y+dy
        case(4)
            call Index1D(k1, n, i, j)
            call Index1D(k2, n, i, j1)
            write(1,fmt) x, interpol_v(k1), x+dx, interpol_v(k2)
        case(5)
            call Index1D(k1, n, i, j1)
            call Index1D(k2, n1, i, j)
            write(1,fmt) x+dx, interpol_v(k1), interpol_h(k2), y
        case(6)
            call Index1D(k1, n, i, j)
            call Index1D(k2, n1, i, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y

            call Index1D(k3, n, i, j1)
            call Index1D(k4, n1, i1, j)
            write(1,fmt) x+dx, interpol_v(k3), interpol_h(k4), y+dy
        case(7)
            call Index1D(k1, n1, i, j)
            call Index1D(k2, n1, i1, j)
            write(1,fmt) interpol_h(k1), y, interpol_h(k2), y+dy
        case(8)
            call Index1D(k1, n, i, j)
            call Index1D(k2, n1, i, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y
        case(9)
            call Index1D(k1, n, i, j)
            call Index1D(k2, n1, i, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y
        case(10)
            call Index1D(k1, n, i, j)
            call Index1D(k2, n, i1, j)
            write(1,fmt) interpol_h(k1), y, interpol_h(k2), y+dy
        case(11)
            call Index1D(k1, n, i, j)
            call Index1D(k2, n1, i1, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y+dy

            call Index1D(k3, n, i, j1)
            call Index1D(k4, n1, i, j)
            write(1,fmt) x+dx, interpol_v(k3), interpol_h(k4), y
        case(12)
            call Index1D(k1, n, i, j1)
            call Index1D(k2, n1, i, j)
            write(1,fmt) x+dx, interpol_v(k1), interpol_h(k2), y
        case(13)
            call Index1D(k1, n, i, j)
            call Index1D(k2, n, i, j1)
            write(1,fmt) x, interpol_v(k1), x+dx, interpol_v(k2)
        case(14)
            call Index1D(k1, n, i, j1)
            call Index1D(k2, n1, i1, j)
            write(1,fmt) x+dx, interpol_v(k1), interpol_h(k2), y+dy
        case(15)
            call Index1D(k1, n, i, j)
            call Index1D(k2, n1, i1, j)
            write(1,fmt) x, interpol_v(k1), interpol_h(k2), y+dy
        
        case default

    end select


end subroutine constructEdges
