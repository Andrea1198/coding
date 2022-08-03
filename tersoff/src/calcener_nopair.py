from numba import jit,config, njit, threading_layer
from numpy import  zeros, array, rint, int32, cos, sin, sqrt, pi, exp

# def randomwalk(x, y, z, N, delta):
#     x  += delta*random.normal(0., 1., N)
#     y  += delta*random.normal(0., 1., N)
#     z  += delta*random.normal(0., 1., N)
# config.THREADING_LAYER = 'threadsafe'

@jit(nopython=True,cache=True, fastmath=True)
def calcener(rx, ry, rz, sp, R, R2, S2, pRSr, C2, H, D2, DR2, X, BETAN, N2R, mu, lam, A, B, Lx, Ly, Lz, mx, my, mz, N, n_) :
    # Zeroing
    fc              = zeros(10)
    fA              = zeros(10)
    fR              = zeros(10)
    # dfA             = zeros(10)
    # dfR             = zeros(10)
    # dfc             = zeros(10)
    dx              = zeros(10)
    dy              = zeros(10)
    dz              = zeros(10)
    r               = zeros(10)
    dr              = zeros((10, 3))
    rr              = zeros(10)
    cosThetaijk     = zeros(10)
    h_cosThetaijk   = zeros(10)
    gThetaijk       = zeros(10)
    spq             = zeros(10, dtype=int32)
    atom            = zeros(10, dtype=int32)
    # fTot            = zeros(3)
    # dgThetaijk      = zeros(10)
    # dCosThetaijk    = zeros((10,3))

    rcx     = zeros(N)
    rcy     = zeros(N)
    rcz     = zeros(N)
    # enep    = 0.
    # virp    = 0.
    # etxx    = zeros(N)
    # stxx    = zeros(N)
    # styy    = zeros(N)
    # stzz    = zeros(N)
    # stxy    = zeros(N)
    # stxz    = zeros(N)
    # styz    = zeros(N)
    e1xx    = zeros(N)
    # s1xx    = zeros(N)
    # s1yy    = zeros(N)
    # s1zz    = zeros(N)
    # s1xy    = zeros(N)
    # s1xz    = zeros(N)
    # s1yz    = zeros(N)
    # fx      = zeros(N)
    # fy      = zeros(N)
    # fz      = zeros(N)

    ncells  = mx*my*mz

    np      = zeros(ncells, dtype=int32)
    indc    = zeros(N, dtype=int32)
    s1p     = zeros(N, dtype=int32)

    for i in range(N):
        vcx=int(mx*(rx[i]+0.5))
        vcy=int(my*(ry[i]+0.5))
        vcz=int(mz*(rz[i]+0.5))
        c = mz*(my*vcx+vcy)+vcz
        indc[i]=c
        np[c] += 1

    indp = zeros(ncells+1, dtype=int32)

    for c in range(0,ncells) :
        indp[c+1] = indp[c] + np[c]

    indcp= zeros(N, dtype=int32)

    for i in range(N):
        c=indc[i]
        rcx[indp[c]] = (rx[i]+0.5)*Lx
        rcy[indp[c]] = (ry[i]+0.5)*Ly
        rcz[indp[c]] = (rz[i]+0.5)*Lz
        s1p[indp[c]] = sp[i]
        indcp[indp[c]] = i
        indp[c] += 1

    # need to reconstruct index list
    indp[:-1] -= np[:] 
    # indp[0]=0
    # for c in range(0,ncells) :
    #     indp[c+1] = indp[c] + np[c]
    # indexing neighbour cells of selected one
    vcx1 = zeros(27, dtype=int32)
    vcy1 = zeros(27, dtype=int32)
    vcz1 = zeros(27, dtype=int32)
    k = 0
    for i in range(-1,2):
        for j in range(-1,2):
            for l in range(-1,2):
                vcx1[k] = i
                vcy1[k] = j
                vcz1[k] = l
                k += 1
    ### PASSAGGI: ###
    # Selezionare solo gli atomi dentro alla sfera S
    # Memorizzare i termini che ricorrono più volte, una lista per ogni atomo i
    enep = 0
    for vcx in range(mx):
        for vcy in range(my):
            for vcz in range(mz):
                c = mz*(my*vcx+vcy)+vcz
                for i in range(indp[c],indp[c+1]):
                    # Coeff i dependent
                    ei      = 0.
                    spi     = s1p[i]
                    c2      = C2[spi]
                    d2      = D2[spi]
                    dr2     = DR2[spi]
                    h       = H[spi]
                    n       = n_[spi]
                    betan   = BETAN[spi]
                    n2r     = N2R[spi]
                    q = 0
                    # loop over particles inside selected cells (central+neighbours)
                    for k in range(27) :
                        wcx=vcx + vcx1[k]
                        wcy=vcy + vcy1[k]
                        wcz=vcz + vcz1[k]
                        # Periodic boundary conditions
                        shiftx = 0.
                        if (wcx == -1) :
                            shiftx =-Lx
                            wcx = mx-1
                        elif (wcx==mx) :
                            shiftx = Lx
                            wcx = 0
                        shifty = 0.
                        if (wcy == -1) :
                            shifty =-Ly
                            wcy = my-1
                        elif (wcy==my) :
                            shifty = Ly
                            wcy = 0
                        shiftz = 0.
                        if (wcz == -1) :
                            shiftz =-Lz
                            wcz = mz-1
                        elif (wcz==mz) :
                            shiftz = Lz
                            wcz = 0
                        c1 = mz*(my*wcx+wcy)+wcz
                        for j in range(indp[c1],indp[c1+1]):
                            dx[q]   = rcx[i]-(rcx[j] + shiftx)
                            dy[q]   = rcy[i]-(rcy[j] + shifty)
                            dz[q]   = rcz[i]-(rcz[j] + shiftz)
                            r2      = dx[q]*dx[q] + dy[q]*dy[q] + dz[q]*dz[q]
                            spq[q]  = s1p[j]
                            index_ij= spq[q]+spi
                            if r2 <= S2[index_ij] and r2 > 0.1:
                                atom[q] = j
                                r[q]    = sqrt(r2)
                                rr[q]   = 1./r[q]
                                dr[q]   = array([dx[q], dy[q], dz[q]])*rr[q]
                                fA[q]   = -B[index_ij]*exp(-mu[index_ij]*r[q])
                                fR[q]   = A[index_ij]*exp(-lam[index_ij]*r[q])
                                # dfA[q]  = -mu[index_ij]*fA[q]
                                # dfR[q]  = -lam[index_ij]*fR[q]
                                if r2 > R2[index_ij]:
                                    a       = pRSr[index_ij]*(r[q]-R[index_ij])
                                    fc[q]   = 0.5 + 0.5*cos(a)
                                    # dfc[q]  = -0.5*sin(a)*pRSr[index_ij]
                                else:
                                    fc[q]   = 1.
                                    # dfc[q]  = 0.
                                q+=1
                        # ond cicle j
                    #end cicle k
                    # Now we have q particles in the S sphere
                    # if q >= 10: return None
                    gThetai = 1. + c2*dr2
                    for next_j1 in range(q):
                        next_j  = next_j1 % q
                        index_ij= spq[next_j] + spi
                        # f3  = 0.
                        # f4  = zeros(3)
                        # f5  = zeros(3)
                        # Calcolo preliminare di zetaij per bij
                        zetaij  = 0.
                        for next_k1 in range(next_j+1, next_j+q):
                            next_k  = next_k1 % q
                            # cos(thetaijk)
                            rij_Scalar_rik          = dx[next_j]*dx[next_k] + dy[next_j]*dy[next_k] + dz[next_j]*dz[next_k]
                            rrij_rrik               = rr[next_j]*rr[next_k]
                            cosThetaijk[next_k]     = rij_Scalar_rik * rrij_rrik                # Memorizzare in un vettore perchè riutilizzata dopo
                            h_cosThetaijk[next_k]   = h - cosThetaijk[next_k]                   # Memorizzare in un vettore perchè riutilizzata dopo
                            gThetaijk_den           = 1./(d2 + h_cosThetaijk[next_k] * h_cosThetaijk[next_k])
                            gThetaijk[next_k]       = gThetai - c2*gThetaijk_den                # Memorizzare in un vettore perchè riutilizzata dopo
                            zetaij                 += fc[next_k]*gThetaijk[next_k]
                            # dgThetaijk[next_k]      = -2.*c2*h_cosThetaijk[next_k] * gThetaijk_den*gThetaijk_den
                            # dCosThetaijk[next_k]    = rr[next_j] * (dr[next_k] - dr[next_j] * cosThetaijk[next_k])

                        bZetaijn    = 1.+betan*(zetaij**n)
                        bij         = X[index_ij]*(bZetaijn)**n2r
                        # dbij        = -0.5*X[index_ij]*(bZetaijn**(n2r-1.))*betan*(zetaij**(n-1.))
                        # for next_k1 in range(next_j+1, next_j+q):
                            # next_k  = next_k1 % q
                            # index_ik= spq[next_k] + spi
                            # cos(thetaijk)
                            # rij_Scalar_rik  = dx[next_j]*dx[next_k] + dy[next_j]*dy[next_k] + dz[next_j]*dz[next_k]
                            # rrij_rrik       = rr[next_j]*rr[next_k]
                            # cosThetaijk     = rij_Scalar_rik * rrij_rrik
                            # h_cosThetaijk   = h-cosThetaijk
                            # gThetaijk_den   = 1./(d2  +  h_cosThetaijk * h_cosThetaijk)
                            # gThetaijk       = gThetai - c2*gThetaijk_den

                            # Calcolo di zetaik
                            # zetaik = 0.
                            # for next_l1 in range(next_k+1, next_k+q):
                                # next_l  = next_l1 % q
                                # rik_Scalar_ril  = dx[next_k]*dx[next_l] + dy[next_k]*dy[next_l] + dz[next_k]*dz[next_l]
                                # rrik_rril       = rr[next_k]*rr[next_l]
                                # cosThetaikl     = rik_Scalar_ril * rrik_rril
                                # h_cosThetaikl   = h-cosThetaikl
                                # gThetaikl_den   = 1./(d2  +  h_cosThetaikl * h_cosThetaikl)
                                # gThetaikl       = gThetai - c2*gThetaikl_den
                                # zetaik         += fc[next_l]*gThetaikl


                            # bZetaikn    = 1.+betan*(zetaik**n)
                            # dbik        = -0.5*X[index_ik]*(bZetaikn**(n2r-1.))*betan*(zetaik**(n-1.))

                            # f34 = fc[next_k] * fA[next_k] * dbik
                            # f3 += dfc[next_j] * gThetaijk[next_k] * f34
                            # f4 += fc[next_j] * dgThetaijk[next_k] * dCosThetaijk[next_k] * f34
                            # f5 += fc[next_k] * fc[next_j] * fA[next_j] * dbij * dgThetaijk[next_k] * dCosThetaijk[next_k]

                        ei     += fc[next_j]*(fR[next_j] + bij*fA[next_j])
                        # f1      = dfc[next_j]*(fR[next_j] + bij*fA[next_j])  # Control - bij or +bij
                        # f2      = fc[next_j]*(dfR[next_j] + bij*dfA[next_j])
                        # fTot    = 0.5*((f1 + f2 + f3)*dr[next_j] + f4 + f5)
                        # fx[indcp[atom[next_j]]] += fTot[0]
                        # fy[indcp[atom[next_j]]] += fTot[1]
                        # fz[indcp[atom[next_j]]] += fTot[2]
                        # fx[indcp[i]] -= fTot[0]
                        # fy[indcp[i]] -= fTot[1]
                        # fz[indcp[i]] -= fTot[2]

                    # Fine ciclo su next_j
                # Fine ciclo sulle celle (k)
                e1xx[indcp[i]] += ei
                mod_r = rcx[indcp[i]]**2 + rcy[indcp[i]]**2 + rcz[indcp[i]]**2
                mod_r = sqrt(mod_r)
                # virp += fx[indcp[i]] * rx[indcp[i]]
                # virp += fy[indcp[i]] * ry[indcp[i]]
                # virp += fz[indcp[i]] * rz[indcp[i]]
                # Fine ciclo su i
            # Fine ciclo su vcz
        # Fine ciclo su vcy
    # Fine ciclo su vcx
    for i in range(N):
        enep += 0.5*e1xx[i]
    # print(enep)
    return enep
