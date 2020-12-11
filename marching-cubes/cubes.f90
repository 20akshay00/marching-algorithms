!reference for tritable; http://paulbourke.net/geometry/polygonise/

program marching_cubes

    use glob
    implicit none
    integer :: i, j, k, l, l1, signcode
    integer, dimension(:), allocatable :: signs !vertex sign wrt. isovalue at each gridpoint
    character(len =10) :: fmt = "(*(1x,g0))"
    character(len = 200) :: lines

    call gridspecInput() !user input to get init data
    call dataGen() !generates function values at gridpoints
    call triTableInit() !initialise tritable values

    allocate(signs(nx*ny*nz))
    numpoints = 0; !counts number of points written in file, for vtk formatting

    do l = 1, nx*ny*nz !calculate signs of each vertex
        if (fgrid(l) >= sigma) then
            signs(l) = 1;
        else 
            signs(l) = 0;
        end if 
    end do

    do i = 1, ny !calculate interpolants of edges along x
        do j = 1, nx-1
            do k = 1, nz
            call Index1D(l, nx, ny, i, j, k)
            call Index1D(l1, nx-1, ny, i, j, k)

            call LinearInterpolation(xgrid(j), fgrid(l), xgrid(j+1), fgrid(l+1), xpol(l1), sigma)
            end do 
        end do
    end do 

    do i = 1, ny-1 !calculate interpolants of edges along y
        do j = 1, nx
            do k = 1, nz
            call Index1D(l, nx, ny, i, j, k)
            call Index1D(l1, nx, ny, i+1, j, k)

            call LinearInterpolation(ygrid(i), fgrid(l), ygrid(i+1), fgrid(l1), ypol(l), sigma)
            end do
        end do
    end do 

    do i = 1, ny !calculate interpolants of edges along z
        do j = 1, nx
            do k = 1, nz-1
            call Index1D(l, nx, ny, i, j, k)
            call Index1D(l1, nx, ny, i, j, k+1)

            call LinearInterpolation(zgrid(k), fgrid(l), zgrid(k+1), fgrid(l1), zpol(l), sigma)
            end do
        end do
    end do 

    open (unit = 1, file = "points.txt", status = "replace")
    do i = 1, ny-1 ! loop to identify cube type and construct vertices
        do j = 1, nx-1
            do k = 1, nz-1
                signcode = 0 !calculate index of cube type; 0-255
                call Index1D(l1, nx, ny, i, j, k)
                signcode = signcode + signs(l1)

                call Index1D(l1, nx, ny, i, j+1, k)
                signcode = signcode + 2*signs(l1)
                
                call Index1D(l1, nx, ny, i, j+1, k+1)
                signcode = signcode + 4*signs(l1)
                
                call Index1D(l1, nx, ny, i, j, k+1)
                signcode = signcode + 8*signs(l1)
                
                call Index1D(l1, nx, ny, i+1, j, k)
                signcode = signcode + 16*signs(l1)

                call Index1D(l1, nx, ny, i+1, j+1, k)
                signcode = signcode + 32*signs(l1)
                
                call Index1D(l1, nx, ny, i+1, j+1, k+1)
                signcode = signcode + 64*signs(l1)
                
                call Index1D(l1, nx, ny, i+1, j, k+1)
                signcode = signcode + 128*signs(l1)

                call constructVertices(signcode, i, j, k, xgrid(j), ygrid(i), zgrid(k))
            end do
        end do
    end do
    close(1)

    !formatting for vtk file 
    open(unit = 1, file = "points.txt", action = "read")
    open(unit = 2, file = "points.vtk", status = "replace")
    
    write(2,fmt) "# vtk DataFile Version 1.0"
    write(2,fmt) "3D triangulation data"
    write(2,fmt) "ASCII"
    write(2,fmt) "DATASET POLYDATA"
    write(2, fmt) "POINTS", numpoints, "float"

    do 
        read(1, '(A)', end = 99) lines
        write(2, '(A)') lines
    end do 

    99 continue 
    close(1, status = 'delete')

    write(2, fmt) "POLYGONS", int(numpoints/3), 4*int(numpoints/3)
    do l = 0, numpoints-1, 3
        write(2, fmt) 3, l, l+1, l+2       
    end do

    close(2)

end program marching_cubes

subroutine LinearInterpolation(x1, f1, x2, f2, x, f) !reverted linear interpolation
    implicit none
    real, intent(in) :: x1, f1, x2, f2, f
    real, intent(out) :: x

    x = x1 + (f-f1)*(x2-x1)/(f2-f1);

end subroutine LinearInterpolation

subroutine Index1D(l, nx, ny, i, j, k) !takes 3D index & #rows, #columns; returns 1D index by flattening row-wise
    implicit none
    integer, intent(in) :: i, j, k, nx, ny
    integer, intent(out) :: l

    l = nx*ny*(k-1) + nx*(i-1) + j 

end subroutine Index1D

subroutine constructVertices(code, i, j, k, x, y, z)
    use glob
    implicit none

    real, intent(in) :: x, y, z
    integer, intent(in) :: code, i, j, k
    character(len =10) :: fmt = "(*(1x,g0))"
    integer, dimension(16) :: edges 
    integer :: l, ind
    ind = 1;
    edges = tritable(code, :)

    do while(edges(ind)/= -1)
        select case(edges(ind)) 
        case(0)
            call Index1D(l, nx-1, ny, i, j, k)
            write(1, fmt) xpol(l), y, z
        case(1)
            call Index1D(l, nx, ny, i, j+1, k)
            write(1, fmt) x+dx, y, zpol(l)
        case(2)
            call Index1D(l, nx-1, ny, i, j, k+1)
            write(1, fmt) xpol(l), y, z+dz
        case(3)
            call Index1D(l, nx, ny, i, j, k)
            write(1, fmt) x, y, zpol(l)
        case(4)
            call Index1D(l, nx-1, ny, i+1, j, k)
            write(1, fmt) xpol(l), y+dy, z
        case(5)
            call Index1D(l, nx, ny, i+1, j+1, k)
            write(1, fmt) x+dx, y+dy, zpol(l)
        case(6)
            call Index1D(l, nx-1, ny, i+1, j, k+1)
            write(1, fmt) xpol(l), y+dy, z+dz
        case(7)
            call Index1D(l, nx, ny, i+1, j, k)
            write(1, fmt) x, y+dy, zpol(l)
        case(8)
            call Index1D(l, nx, ny, i, j, k)
            write(1, fmt) x, ypol(l), z
        case(9)
            call Index1D(l, nx, ny, i, j+1, k)
            write(1, fmt) x+dx, ypol(l), z
        case(10)
            call Index1D(l, nx, ny, i, j+1, k+1)
            write(1, fmt) x+dx, ypol(l), z+dz
        case(11)
            call Index1D(l, nx, ny, i, j, k+1)
            write(1, fmt) x, ypol(l), z+dz
        case default

        end select

        numpoints = numpoints + 1
        ind = ind + 1
    end do

end subroutine constructVertices
