MODULE qties_module
    USE kinds
    !
    IMPLICIT NONE
    !
    INTEGER  :: ndim=0
    INTEGER  :: npart=0
    !
    REAL(dbl)   :: pos(:,:), vel(:,:), acc(:,:)
    REAL(dbl)   :: mass(:)
    !
    LOGICAL :: alloc=.false.
    !
    !-------------------------!
    CONTAINS 
    !-------------------------!
        !
        SUBROUTINE allocate(mass_, npp)
            USE kinds
            !
            IMPLICIT NONE
            INTEGER :: ierr
            REAL(dbl) :: mass_(npp)
            !
            IF (ndim .le. 0 .or. npart .le. 0) CALL error(subname, "ndim or npart not positive integer", ndim + npart)
            ALLOCATE (pos(ndim,npart), vel(ndim,npart), acc(ndim,npart))
            ALLOCATE (mass(npart))
            !
            IF (npp .ne. npart) CALL error(subname, "npart different from npp", npart - ndim)
            !
            mass = mass_
            alloc = .true.
            !
        END SUBROUTINE allocate
        !
        SUBROUTINE deallocate()
            IMPLICIT NONE
            !
            INTEGER :: ierr
            !
            DEALLOCATE (pos, vel, acc, mass)
            alloc=.false.
        END SUBROUTINE deallocate
END MODULE qties_module

