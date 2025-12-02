module defvar
    implicit none
    private
    public :: n1, n2, n3, ncenter, orgx, orgy, orgz, &
              v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, &
              detmat, readcube, calc_dvol

    integer :: n1, n2, n3, ncenter
    real(kind=8) :: orgx, orgy, orgz
    real(kind=8) :: v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z

contains
    real(kind=8) function detmat(mat)
        real(kind=8), intent(in) :: mat(3,3)
        detmat = mat(1,1)*(mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
               - mat(1,2)*(mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
               + mat(1,3)*(mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))
    end function detmat

    subroutine readcube(cubname, cubmat_array)
        character(len=*), intent(in) :: cubname
        real(kind=8), allocatable, intent(out) :: cubmat_array(:,:,:)
        integer :: i, j, ierr

        open(10, file=cubname, status="old", iostat=ierr)
        if (ierr /= 0) error stop "File open failed"

        read(10,*) 
        read(10,*) 
        read(10,*) ncenter, orgx, orgy, orgz
        read(10,*) n1, v1x, v1y, v1z
        read(10,*) n2, v2x, v2y, v2z
        read(10,*) n3, v3x, v3y, v3z

        allocate(cubmat_array(n1,n2,n3), stat=ierr)
        if (ierr /= 0) error stop "Allocation failed"

        do i = 1, ncenter
            read(10,*) 
        end do

        do i = 1, n1
            do j = 1, n2
                read(10,*) cubmat_array(i,j,:)
            end do
        end do
        close(10)
    end subroutine readcube

    subroutine calc_dvol(dvol)
        real(kind=8), intent(out) :: dvol
        real(kind=8) :: mat(3,3)

        mat(:,1) = [v1x, v1y, v1z]
        mat(:,2) = [v2x, v2y, v2z]
        mat(:,3) = [v3x, v3y, v3z]
        dvol = abs(detmat(mat))
    end subroutine calc_dvol
end module defvar

program coulomb_attraction
    use defvar
    use omp_lib
    implicit none
    character(len=200) :: homofile, lumofile
    real(kind=8) :: coulene, dvol, au2eV
    real(kind=8), allocatable :: homomat(:,:,:), lumomat(:,:,:)
    real(kind=8), allocatable :: xcoord(:), ycoord(:), zcoord(:)
    integer :: i, j, k, ii, jj, kk, nthreads, ierr
    real(kind=8) :: t0, t1, r
    logical :: file_exists

    nthreads = omp_get_max_threads()
    call omp_set_num_threads(nthreads)

    au2eV = 27.2114d0

    homofile = "homo.cube"
    inquire(file=homofile, exist=file_exists)
    if (.not. file_exists) then
        write(*,*) "homo.cube not found. Enter the HOMO cube file name:"
        read(*,*) homofile
    end if
    call readcube(homofile, homomat)

    lumofile = "lumo.cube"
    inquire(file=lumofile, exist=file_exists)
    if (.not. file_exists) then
        write(*,*) "lumo.cube not found. Enter the LUMO cube file name:"
        read(*,*) lumofile
    end if
    call readcube(lumofile, lumomat)

    allocate(xcoord(n1), ycoord(n2), zcoord(n3), stat=ierr)
    if (ierr /= 0) error stop "Allocation failed"

    do i = 1, n1
        xcoord(i) = orgx + (i-1)*v1x + (i-1)*v2x + (i-1)*v3x
    end do
    do j = 1, n2
        ycoord(j) = orgy + (j-1)*v1y + (j-1)*v2y + (j-1)*v3y
    end do
    do k = 1, n3
        zcoord(k) = orgz + (k-1)*v1z + (k-1)*v2z + (k-1)*v3z
    end do

    coulene = 0.0d0
    t0 = omp_get_wtime()

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,ii,jj,kk,r) REDUCTION(+:coulene)
    !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                if (abs(homomat(i,j,k)) < 1D-6) cycle

                do kk = 1, n3
                    do jj = 1, n2
                        do ii = 1, n1
                            if (i /= ii .or. j /= jj .or. k /= kk) then
                                r = sqrt((xcoord(i)-xcoord(ii))**2 + &
                                        (ycoord(j)-ycoord(jj))**2 + &
                                        (zcoord(k)-zcoord(kk))**2)
                                coulene = coulene + (homomat(i,j,k)**2) * (lumomat(ii,jj,kk)**2) / r
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    t1 = omp_get_wtime()

    call calc_dvol(dvol)
    coulene = coulene * dvol**2

    write(*,"(' Coulomb attractive energy:',f12.6,' a.u.  (',f12.6,' eV )')") coulene, coulene*au2eV
    write(*,"(' Calculation time:',f8.2,' seconds')") t1-t0

    deallocate(homomat, lumomat, xcoord, ycoord, zcoord)
end program coulomb_attraction