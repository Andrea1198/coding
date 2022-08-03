SUBROUTINE update(x, v, a)
    use kinds, only : dbl
    use qties, only : x, v, a, dt, ndim
    implicit none
    !
    x = x + v*dt
    v = v + a*dt
END SUBROUTINE
