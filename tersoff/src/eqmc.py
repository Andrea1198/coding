from numpy import sum, zeros, rint, sqrt, random, int32, exp
from src.calcener_MC import calcener_MC
from src.calcener_nopair import calcener
#import tqdm
def eqmc(self, N, kt, nstep, dt, freq, mode, dir, delta, dV0, pres):
    # initializing counters and constants
    dth=0.5*dt
    fout=open(dir+'equi.txt','a')
    t    = 0
    self.ept  = 0.
    self.ekt  = 0.
    self.pres = 0.
    rx = zeros(N)
    ry = zeros(N)
    rz = zeros(N)
    dacc = 0
    moved= 0
    vacc = 0
    movev= 0

    fout.write("# starting eqmc trajectory  Nb = %7d\n" % sum(self.sp) )
    fout.write( "    'pas'    'enep'    'enek'      'enet'      'vcm'\n")
    # initial evaluation of forces, energies and virials
    # rx  =  self.x/self.Lx
    # rx -= rint(rx)
    # ry  =  self.y/self.Ly
    # ry -= rint(ry)
    # rz  =  self.z/self.Lz
    # rz -= rint(rz)
    # (enep, virial, self.etxx, self.stxx, self.styy, self.stzz, self.stxy, self.stxz, self.styz, self.fx, self.fy, self.fz) = \
    # self.fx, self.fy, self.fz, enep = calcener( rx, ry, rz, self.sp, self.R, self.R2, self.S2, self.pRSr, self.c2, self.h, self.d2, self.dr2, self.X, self.betan, self.n2r, self.mu, self.lam, self.A, \
    # self.B, self.Lx, self.Ly, self.Lz, self.mx, self.my, self.mz, N, self.n)
    # fout.write (" %8.3f %9.4f \n" % (0., enep/N ) )


    for pas in range(nstep) :
        # random particle position
        mover = random.randint(0,N+2)
        rx[0:N] = self.x[0:N] / self.Lx
        rx[0:N]-= rint(rx[0:N])
        ry[0:N] = self.y[0:N] / self.Ly
        ry[0:N]-= rint(ry[0:N])
        rz[0:N] = self.z[0:N] / self.Lz
        rz[0:N]-= rint(rz[0:N])
        if mover < N:
            moved += 1
            # mover = random.randint(0,N,n_movers)
            enep_in = calcener_MC(rx, ry, rz, self.sp, self.R, self.R2, self.S2, self.pRSr, self.c2, self.h, self.d2, self.dr2, self.X, self.betan, self.n2r, self.mu, self.lam, self.A, \
                self.B, self.Lx, self.Ly, self.Lz, self.mx, self.my, self.mz, N, self.n, mover)
            dx = delta*random.rand()*2 - delta
            rx[mover] += dx
            dy = delta*random.rand()*2 - delta
            ry[mover] += dy
            dz = delta*random.rand()*2 - delta
            rz[mover] += dz

            # Periodic boundary conditions
            if rx[mover] > 0.5:
                rx[mover] -= 1
            elif rx[mover] < -0.5:
                rx[mover] += 1
            if ry[mover] > 0.5:
                ry[mover] -= 1
            elif ry[mover] < -0.5:
                ry[mover] += 1
            if rz[mover] > 0.5:
                rz[mover] -= 1
            elif rz[mover] < -0.5:
                rz[mover] += 1
            # rx[rx > 0.5] -= 1
            # rx[rx < -0.5] += 1

            # ry[ry > 0.5] -= 1
            # ry[ry < -0.5] += 1

            # rz[rz > 0.5] -= 1
            # rz[rz < -0.5] += 1

            enep_fin = calcener_MC(rx, ry, rz, self.sp, self.R, self.R2, self.S2, self.pRSr, self.c2, self.h, self.d2, self.dr2, self.X, self.betan, self.n2r, self.mu, self.lam, self.A, \
                self.B, self.Lx, self.Ly, self.Lz, self.mx, self.my, self.mz, N, self.n, mover)
            enep_diff = enep_in - enep_fin
            accep = exp(enep_diff/kt)
            if random.rand() < accep:
                self.x[mover] += dx * self.Lx
                self.y[mover] += dy * self.Ly
                self.z[mover] += dz * self.Lz
                dacc += 1
        # Change volume
        else:
            movev += 1
            # rx[0:N] = self.x[0:N] / self.Lx
            # rx[0:N]-= rint(rx[0:N])
            # ry[0:N] = self.y[0:N] / self.Ly
            # ry[0:N]-= rint(ry[0:N])
            # rz[0:N] = self.z[0:N] / self.Lz
            # rz[0:N]-= rint(rz[0:N])
            enep_in = calcener(rx, ry, rz, self.sp, self.R, self.R2, self.S2, self.pRSr, self.c2, self.h, self.d2, self.dr2, self.X, self.betan, self.n2r, self.mu, self.lam, self.A, \
                self.B, self.Lx, self.Ly, self.Lz, self.mx, self.my, self.mz, N, self.n)\

            temp = random.rand()*2*dV0 - dV0

            # Change volume uniformely
            coef = pow(1. + temp, 1/3) 
            Lx = self.Lx * coef
            Ly = self.Ly * coef
            Lz = self.Lz * coef
            # Change only one direction
            # direc = random.randint(2)
            # if direc == 0:
            #     self.Lx += temp * self.volume / (self.Ly*self.Lz)
            # elif direc == 1:
            #     self.Ly += temp * self.volume / (self.Lx*self.Lz)
            # else: # direc == 2
            #     self.Lz += temp * self.volume / (self.Ly*self.Lx)
            enep_fin = calcener(rx, ry, rz, self.sp, self.R, self.R2, self.S2, self.pRSr, self.c2, self.h, self.d2, self.dr2, self.X, self.betan, self.n2r, self.mu, self.lam, self.A, \
                self.B, Lx, Ly, Lz, self.mx, self.my, self.mz, N, self.n)
            if enep_fin != None:
                enep_diff = enep_in-enep_fin
                accep = coef ** N * exp((enep_diff - pres*dV0)/kt)
                if random.rand() < accep:
                    self.Lx = Lx
                    self.Ly = Ly
                    self.Lz = Lz 
                    self.volume = Lx * Ly * Lz
                    vacc += 1
                
        if pas%(freq) == 0 and mode == 1:
            self.calcgdr(N, rx, ry, rz)

        # if pas%(freq*100) == 0 and mode == 2:
        #     # andersen thermostats: velocity rescaling
        #     pstd=sqrt(self.m*kt)
        #     self.px[0:N] = pstd[0:N]*random.normal(0., 1., N)
        #     self.py[0:N] = pstd[0:N]*random.normal(0., 1., N)
        #     self.pz[0:N] = pstd[0:N]*random.normal(0., 1., N)
        #     mtot  = sum(self.m)
        #     vcmx  = sum(self.px)/mtot
        #     vcmy  = sum(self.py)/mtot
        #     vcmz  = sum(self.pz)/mtot
        #     self.px[0:N] -= self.m[0:N]*vcmx
        #     self.py[0:N] -= self.m[0:N]*vcmy
        #     self.pz[0:N] -= self.m[0:N]*vcmz
        #     fout.write("# velocities sampled from maxwell distribution at timestep %d\n" % pas)
            
        t   += dt
        enek = 0.5*sum ( (self.px[0:N]**2 + self.py[0:N]**2 + self.pz[0:N]**2)/self.m[0:N] )
        self.ekt += enek
        self.ept += enep_in
        # self.pres+= virial
        if (pas+1)%(freq*10)==0 :
            fout.write (" %8.3f %9.4f \n" % \
            (t, enep_in/N) )
    # end of md run
    g=3*N-3
    fout.write( "# ending eqmc trajectory  <ep>=%6.1f  <ek>=%7.4f   T=%8.3f  P=%8.3f\n" % ( self.ept/nstep, self.ekt/nstep, 2.*self.ekt/(nstep*g), self.pres/nstep) )
    fout.write( "# acceptance rates displacements=%6.2f  volume=%6.2f   ratio D/VT=%6.2f\n" % ( dacc/moved, vacc/movev, moved/movev) )
    fout.close()
