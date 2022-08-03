SUBROUTINE error(subname, error_string, ierr)
    USE constants,  ONLY : stdout
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)         :: ierr
    CHARACTER(256), INTENT(IN)  :: error_string
    CHARACTER(32), INTENT(IN)   :: subname
    !
    WRITE(stdout,*) "######################################################################"
    WRITE(stdout,*) "Program stopped in ", subname, " because of ", error_string
    WRITE(stdout,*) "Error code : ", ierr
    WRITE(stdout,*) "######################################################################"
    stop
END SUBROUTINE error
