MODULE qties
    INTEGER                 :: ndim=0
    REAL(dbl), ALLOCATABLE  :: x(:), v(:), a(:)
    LOGICAL     :: alloc=.false.

    PUBLIC :: ndim
    PUBLIC :: x, v, a 
    PUBLIC :: alloc
    PUBLIC :: system_allocate
    PUBLIC :: system_deallocate

CONTAINS
    SUBROUTINE system_allocate()
        implicit none
        integer :: ierr
        character(15) :: "system_allocate"
        !
        if (alloc) call errore(subname, "system already allocated", 1)
        ndim=3
        allocate(x(ndim), v(ndim), a(ndim), stat=ierr)
        if (ierr /= 0) call errore(subname, "error in allocating x, v, a", ierr)
    end subroutine system_allocate

    subroutine system_deallocate()
        implicit none
        character(17) :: subname="system_deallocate"
        !
        if ( .not. alloc ) call errore(subname, "system already deallocated", 2)
        ndim = 0
        deallocate(x, v, a)
    end subroutine system_deallocate
