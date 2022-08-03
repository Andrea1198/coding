from numba import jit,config, njit, threading_layer
from numpy import  zeros, array, rint, int32, cos, sin, sqrt, pi, exp

# def randomwalk(x, y, z, N, delta):
#     x  += delta*random.normal(0., 1., N)
#     y  += delta*random.normal(0., 1., N)
#     z  += delta*random.normal(0., 1., N)
# config.THREADING_LAYER = 'threadsafe'

# @jit(nopython=True,cache=True, fastmath=True)
def calcener_MC(rx, ry, rz, sp, R, R2, S2, pRSr, C2, H, D2, DR2, X, BETAN, N2R, mu, lam, A, B, Lx, Ly, Lz, mx, my, mz, N, n_, mover) :
    # Zeroing
    fc              = zeros(10)
    fA              = zeros(10)
    fR              = zeros(10)
    dx              = zeros(10)
    dy              = zeros(10)
    dz              = zeros(10)
    r               = zeros(10)
    rr              = zeros(10)
    spq             = zeros(10, dtype=int32)
    atom            = zeros(10, dtype=int32)

    rcx     = zeros(N)
    rcy     = zeros(N)
    rcz     = zeros(N)

    ncells  = mx*my*mz

    np      = zeros(ncells, dtype=int32)
    indc    = zeros(N, dtype=int32)
    s1p     = zeros(N, dtype=int32)

    # for i in range(N):
    #     vcx=int(mx*(rx[i]+0.5))
    #     vcy=int(my*(ry[i]+0.5))
    #     vcz=int(mz*(rz[i]+0.5))
    #     c = mz*(my*vcx+vcy)+vcz
    #     indc[i]=c
    #     np[c] += 1

    # indp = zeros(ncells+1, dtype=int32)

    # for c in range(0,ncells) :
    #     indp[c+1] = indp[c] + np[c]

    # indcp= zeros(N, dtype=int32)

    rcx = (rx+0.5)*Lx
    rcy = (ry+0.5)*Ly
    rcz = (rz+0.5)*Lz

    # # need to reconstruct index list
    # indp[:-1] -= np[:] 
    # # indp[0]=0
    # # for c in range(0,ncells) :
    # #     indp[c+1] = indp[c] + np[c]
    # # indexing neighbour cells of selected one
    # vcx1 = zeros(27, dtype=int32)
    # vcy1 = zeros(27, dtype=int32)
    # vcz1 = zeros(27, dtype=int32)
    # k = 0
    # for i in range(-1,2):
    #     for j in range(-1,2):
    #         for l in range(-1,2):
    #             vcx1[k] = i
    #             vcy1[k] = j
    #             vcz1[k] = l
    #             k += 1
    # ### PASSAGGI: ###
    # # Selezionare solo gli atomi dentro alla sfera S
    # # Memorizzare i termini che ricorrono più volte, una lista per ogni atomo i
    # vcx=int(mx*(rx[mover]+0.5))
    # vcy=int(my*(ry[mover]+0.5))
    # vcz=int(mz*(rz[mover]+0.5))
    ei      = 0.

    spi     = sp[mover]
    c2      = C2[spi]
    d2      = D2[spi]
    dr2     = DR2[spi]
    h       = H[spi]
    n       = n_[spi]
    betan   = BETAN[spi]
    n2r     = N2R[spi]
    # q = 0
    # # loop over particles inside selected cells (central+neighbours)
    # for k in range(27):
    #     wcx=vcx + vcx1[k]
    #     wcy=vcy + vcy1[k]
    #     wcz=vcz + vcz1[k]
    #     # Periodic boundary conditions
    #     shiftx = 0.
    #     if (wcx == -1) :
    #         shiftx =-Lx
    #         wcx = mx-1
    #     elif (wcx==mx) :
    #         shiftx = Lx
    #         wcx = 0
    #     shifty = 0.
    #     if (wcy == -1) :
    #         shifty =-Ly
    #         wcy = my-1
    #     elif (wcy==my) :
    #         shifty = Ly
    #         wcy = 0
    #     shiftz = 0.
    #     if (wcz == -1) :
    #         shiftz =-Lz
    #         wcz = mz-1
    #     elif (wcz==mz) :
    #         shiftz = Lz
    #         wcz = 0
    #     c1 = mz*(my*wcx+wcy)+wcz
    q = 0
    for j in range(N):
        dx[q]   = rcx[mover]-(rcx[j])
        dy[q]   = rcy[mover]-(rcy[j])
        dz[q]   = rcz[mover]-(rcz[j])

        if dx[q] > Lx/2. :
            dx[q] -= Lx
        elif dx[q] < -Lx/2. :
            dx[q] += Lx
        if dy[q] > Ly/2. :
            dy[q] -= Ly
        elif dy[q] < -Ly/2. :
            dy[q] += Ly
        if dz[q] > Lz/2. :
            dz[q] -= Lz
        elif dz[q] < -Lz/2. :
            dz[q] += Lz

        r2      = dx[q]*dx[q] + dy[q]*dy[q] + dz[q]*dz[q]
        spq[q]  = sp[j]
        index_ij= spq[q]+spi
        if r2 <= S2[index_ij] and r2 > 0.1:
            atom[q] = j
            r[q]    = sqrt(r2)
            rr[q]   = 1./r[q]
            fA[q]   = -B[index_ij]*exp(-mu[index_ij]*r[q])
            fR[q]   = A[index_ij]*exp(-lam[index_ij]*r[q])
            if r2 > R2[index_ij]:
                a       = pRSr[index_ij]*(r[q]-R[index_ij])
                fc[q]   = 0.5 + 0.5*cos(a)
            else:
                fc[q]   = 1.
            q+=1
        # ond cicle j
    #end cicle k
    if q >= 10: return None
    # Now we have q particles in the S sphere
    gThetai = 1. + c2*dr2
    for next_j1 in range(q):
        next_j  = next_j1 % q
        index_ij= spq[next_j] + spi
        # Calcolo preliminare di zetaij per bij
        zetaij  = 0.
        for next_k1 in range(next_j+1, next_j+q):
            next_k  = next_k1 % q
            # cos(thetaijk)
            rij_Scalar_rik          = dx[next_j]*dx[next_k] + dy[next_j]*dy[next_k] + dz[next_j]*dz[next_k]
            rrij_rrik               = rr[next_j]*rr[next_k]
            cosThetaijk             = rij_Scalar_rik * rrij_rrik
            h_cosThetaijk           = h - cosThetaijk
            gThetaijk_den           = 1./(d2 + h_cosThetaijk * h_cosThetaijk)
            gThetaijk               = gThetai - c2*gThetaijk_den
            zetaij                 += fc[next_k]*gThetaijk

        bZetaijn    = 1.+betan*(zetaij**n)
        bij         = X[index_ij]*(bZetaijn)**n2r
        ei     += fc[next_j]*(fR[next_j] + bij*fA[next_j])

        # Fine ciclo su next_j
    # Fine ciclo sulle celle (k)
    # mod_r = rcx[indcp[mover]]**2 + rcy[indcp[mover]]**2 + rcz[indcp[mover]]**2
    # mod_r = sqrt(mod_r)
                # Fine ciclo su i
            # Fine ciclo su vcz
        # Fine ciclo su vcy
    # Fine ciclo su vcx
    # print(enep)
    return ei

