program marching_squares
    implicit none
    use marching_variables

    call gridspecInput(xmin, xmax, dx, m, ymin, ymax, dy, n, xgrid, ygrid, fgrid)
    call dataGen(xgrid, m, ygrid, n, fgrid)

end program marching_squares

subroutine LinearInterpolation(x1, f1, x2, f2, x, f)
    implicit none
    real, intent(in) :: x1, f1, x2, f2, f
    real, intent(out) :: x

    x = x1 + (f-f1)*(x2-x1)/(f2-f1)
end subroutine LinearInterpolation

subroutine Index1D(i, j, m, k) !takes 2D index & #rows; returns 1D index by flattening row-wise
    implicit none
    integer, intent(in) :: i, j, m
    integer, intent(out) :: k

    k = m*(i-1) + j

end subroutine Index1D

subroutine Index2D(k, m, i, j) !reverse flatten, returns 2D index using 1D index 
    implicit none
    integer, intent(in) :: k, m
    integer, intent(out) :: i, j
    
    i = 1 + k/m 
    j = modulo(k, m)

end subroutine Index2D
