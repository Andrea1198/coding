class TER :

    def __init__(self, l1, l2, l3):
        from numpy import zeros, ones, sqrt, int32, float32, loadtxt, array, pi
        #-start-------------------------
        #TERSOFF
        self.ma     = 16.    # O : 15.9999 amu
        self.mb     = 28.    # Si : 28.0855
        file        = './files/input/coeff_SiO_tersoff.txt'
        n_atoms = loadtxt(file, delimiter=',', max_rows=1, dtype=int32)
        combine = 3
        self.A      = zeros(combine)
        self.B      = zeros(combine)
        self.lam    = zeros(combine)
        self.mu     = zeros(combine)
        self.beta   = zeros(n_atoms)
        self.c      = zeros(n_atoms)
        self.d      = zeros(n_atoms)
        self.h      = zeros(n_atoms)
        self.X      = zeros(combine)
        self.S      = zeros(combine)
        self.R      = zeros(combine)
        self.n      = zeros(n_atoms)
        self.A[::2], self.B[::2], self.lam[::2], self.mu[::2], self.beta[:], self.n[:], self.c[:], self.d[:], self.h[:], self.R[::n_atoms], self.S[::n_atoms], \
        x = loadtxt(file, delimiter=',', skiprows=2, max_rows=n_atoms, unpack=True)
        x = x[0]

        self.d2     = self.d**2
        self.c2     = self.c**2
        self.dr2    = 1./self.d2
        self.A[1]   = sqrt(self.A[0]*self.A[2])
        self.B[1]   = sqrt(self.B[0]*self.B[2])
        self.lam[1] = (self.lam[0]+self.lam[2])/2.
        self.mu[1]  = (self.mu[0]+self.mu[2])/2.
        self.R[1]   = sqrt(self.R[0]*self.R[2])
        self.S[1]   = sqrt(self.S[0]*self.S[2])
        self.X      = array([1., x, 1.])

        self.R2     = self.R*self.R
        self.S2     = self.S*self.S
        self.pRSr   = pi/(self.S-self.R)
        self.betan  = self.beta**self.n
        self.n2r    = -1./(2.*self.n)

        # box need be an integer number of rcut
        # Vectors of the Bravais cell needed to generate the crystal

        # Read input file
        self.read_vec_file(filename='./files/input/alphaquartz.csv')

        self.Lx = l1 * self.ax
        self.Ly = l2 * self.by
        self.Lz = l3 * self.cz
        self.rcut   = max(self.S)
        self.mx = int(self.Lx // self.rcut)
        self.my = int(self.Ly // self.rcut)
        self.mz = int(self.Lz // self.rcut)
        self.kb   = 1.             # Boltzmann constant in these units
        # print(self.mx, self.my, self.mz, '\n', self.Lx, self.Ly, self.Lz, '\n', self.rcut)
        # self.lb = self.mz*3
        # self.mx = l1
        # self.my = l2
        # self.mz = l3

        self.npart  = self.nvec*l1*l2*l3
        self.volume = self.Lx*self.Ly*self.Lz
        self.rho    = self.npart/self.volume
        ndim = self.npart
        # N initial guess of average number of particles for dimensioning
        self.x       = zeros( ndim )
        self.y       = zeros( ndim )
        self.z       = zeros( ndim )
        self.m       = zeros( ndim )
        self.sp      = ones( ndim, dtype=int32 )
        self.rx      = zeros( ndim )
        self.ry      = zeros( ndim )
        self.rz      = zeros( ndim )
        self.px      = zeros( ndim )
        self.py      = zeros( ndim )
        self.pz      = zeros( ndim )
        self.fx      = zeros( ndim )
        self.fy      = zeros( ndim )
        self.fz      = zeros( ndim )
        self.etxx    = zeros( ndim )
        self.stxx    = zeros( ndim )
        self.styy    = zeros( ndim )
        self.stzz    = zeros( ndim )
        self.stxy    = zeros( ndim )
        self.stxz    = zeros( ndim )
        self.styz    = zeros( ndim )
        #
        self.kg      = 1024
        self.gcount  = zeros( (self.kg,3) )
        self.integral = zeros((self.kg,4))
        self.ekin    = 0.0
        self.ene     = 0.0
        self.etot    = 0.0
        #
        self.tt      = 0.0
        self.ekt     = 0.0
        self.ept     = 0.0
        self.pres    = 0.0
        #
        rmax = min( (self.Lx, self.Ly, self.Lz) )/2.
        # rmax = self.rcut # to compute gdr in calcener

        # assuming to look up to Lx/2.
        self.r2max = rmax * rmax
        self.ldel = rmax/self.kg

        # assuming 2 different atoms --> 3 pairs

    def calcgdr(self, N, rx, ry, rz):
        from numpy import sqrt, rint, zeros, int32, pi
        for k in range(N-1) :
            j=k+1
            #for j in range(k+1,N) :
            dx = rx[k]-rx[j:N]
            dy = ry[k]-ry[j:N]
            dz = rz[k]-rz[j:N]
            dx-= rint(dx)
            dy-= rint(dy)
            dz-= rint(dz)
            dx = dx*self.Lx
            dy = dy*self.Ly
            dz = dz*self.Lz
            irdf = zeros(N-j, dtype=int32)
            irdf[:] = self.sp[j:N]+self.sp[k]
            r2 = dx*dx + dy*dy + dz*dz
            # using the mask array "b" for speedup
            b = r2 < self.r2max
            lm  = sqrt(r2[b])
            ind = irdf[b]
            #if lm<self.kg :
            for j in range(len(lm)) :
                self.gcount[int(lm[j]/self.ldel),ind[j]]+=2.


    def integrate_gdr(self, g):
        from numpy import zeros, pi, array
        # integral = [O-O,O-Si,Si-Si,Si-O]
        r = self.ldel
        self.integral = zeros(( self.kg, 4 ))
        self.integral[0,0:3] = g[0,:]*r*r*r*0.5
        for i in range(1,self.kg-1) :
            r += self.ldel
            self.integral[i, 0:3] = self.integral[i-1, 0:3] + ((g[i-1, 0:3] + g[i, 0:3])*0.5)*r*r*self.ldel
        self.integral[self.kg-1, 0:3] += g[self.kg-1, :]
        self.integral  *= pi*4.
        self.integral[:,3] = self.integral[:,1]
        rho = array([(self.Na-1)/self.volume, self.Nb/self.volume, (self.Nb-1)/self.volume, self.Na/self.volume])
        self.integral *= rho

    def read_vec_file(self, filename='alphaquartz.csv'):
        from numpy import zeros, loadtxt, int32
        self.nvec       = loadtxt(filename, delimiter=',', max_rows=1, dtype=int32)
        self.xvector    = zeros(self.nvec)
        self.yvector    = zeros(self.nvec)
        self.zvector    = zeros(self.nvec)
        self.spv        = zeros(self.nvec, int32)
        self.ax,self.by,self.cz = loadtxt(filename, delimiter=',', skiprows=1, max_rows=1)
        self.spv[:],self.xvector[:],self.yvector[:],self.zvector[:]= loadtxt(filename, delimiter=',', skiprows=2, unpack=True)
        self.spv = int32(self.spv)

    def writexyz(self, N, filename='conf.xyz'):
        from numpy import savetxt, column_stack, empty, str_, rint
        dx = self.x/self.Lx
        dx-= rint(dx)
        self.x = dx*self.Lx
        dy = self.y/self.Ly
        dy-= rint(dy)
        self.y = dy*self.Ly
        dz = self.z/self.Lz
        dz-= rint(dz)
        self.z = dz*self.Lz
        ar = empty(N,(str_,3))
        sig=1 # in Angstroem for argon
        for i in range(N):
            if(self.sp[i]==0):
                ar[i] = "O "
            else :
                ar[i] = "Si"
        dx *= sig*self.Lx
        dy *= sig*self.Ly
        dz *= sig*self.Lz
        # vx  = self.px/self.m
        # vy  = self.py/self.m
        # vz  = self.pz/self.m
        rout = column_stack( (ar, dx, dy, dz) )
        savetxt(filename, rout, fmt=(' %s',' %s',' %s',' %s'), \
        header=(' %d' % N ), comments=' ' )



    def read_binput(self, N, conf_in='conf_in.b'): #, mom_in='mom_in'):
        #from numpy import loadtxt, float64
        import pickle
        with open(conf_in, 'rb') as ftrj:
            (Nr, Nb, pas) = pickle.load(ftrj)
            if N!=Nr :
                print(' reading %d , %d particle configuration from step %d' % (Nr,Nb, pas) )
            #else :
                print(' ??? reading %d particle configuration expected %d' % (Nr,N) )
            ( self.x,  self.y,  self.z , self.sp) = pickle.load( ftrj)
            ( self.px, self.py, self.pz, self.m ) = pickle.load( ftrj)
        return pas

    def read_input(self, N, conf_in='conf_in.csv'): #, mom_in='mom_in'):
        from numpy import loadtxt, float64, int32
        print('# reading %d particle configuration' % N )
        (self.x, self.y,  self.z , self.sp, self.px, self.py, self.pz, self.m) = loadtxt( conf_in, skiprows=1, unpack=True, delimiter=',')
        self.sp = int32(self.sp)

    def write_binput(self, N, pas, conf_out='conf_in.b'): #, mom_out='mom_in'):
        from numpy import sum
        import pickle
        with open(conf_out, 'wb') as ftrj:
            pickle.dump( (N, sum(self.sp), pas) , ftrj, pickle.HIGHEST_PROTOCOL)
            pickle.dump( ( self.x,  self.y,  self.z , self.sp), ftrj, pickle.HIGHEST_PROTOCOL)
            pickle.dump( ( self.px, self.py, self.pz, self.m ), ftrj, pickle.HIGHEST_PROTOCOL)

    def write_input(self, N, pas, conf_out='conf_out.csv'):
        from numpy import sum, zeros, pi, savetxt, column_stack
        allout = column_stack( (self.x,  self.y,  self.z , self.sp, self.px, self.py, self.pz, self.m) )
        savetxt(conf_out, allout, delimiter=',', header="x,y,z,sp,vx,vy,vz,m" )


    def write_out(self, N, tstep, gdr_out='gdr.out'):
        from numpy import zeros, pi, savetxt, column_stack, sum
        V = zeros(self.kg)
        r = zeros(self.kg)
        g = zeros( (self.kg,3) )
        Nb = sum(self.sp)
        Na = N - Nb
        for lm in range(self.kg) :
            V[lm] = 4./3.*pi*(self.ldel**3)*(3*lm*lm +3*lm + 1);
            r[lm] = (lm+0.5)*self.ldel
        g[:,0] = self.gcount[:,0]/(V*(Na-1)*tstep*self.rho*Na/N);
        g[:,1] = self.gcount[:,1]/(2.*V*Na*Nb*tstep*self.rho/N);
        g[:,2] = self.gcount[:,2]/(V*(Nb-1)*tstep*self.rho*Nb/N);
        self.integrate_gdr(g)
        gout = column_stack( (r, g) )
        savetxt(gdr_out, gout , fmt=('%10.5g ','%12.7g','%12.7g','%12.7g'), \
        header="    'r'     'gaa'     'gab'     'gbb' " )


    def fcc(self):
        from numpy import  random, arange
        # Unitary cell lengths
        ax = self.Lx/self.l1
        ay = self.Ly/self.l2
        az = self.Lz/self.l3
        print( "# lattice parameters a =%.4f  %.4f  %.4f" % (ax, ay, az) )
        print( "# mx = %d   my = %.d  mz = %.d" % (self.l1,self.l2, self.l3) )
        m = self.l1*self.l2*self.l3
        natom = 4*m
        print( "# number of lattice cells m = %d" % m )
        print( "# number of particles  %d " % natom )
        print( "# sides of md-box L = [ %.4f %.4f %.4f ]" % (self.Lx, self.Ly, self.Lz) )
        j  = 0
        xi = 0.
        yi = 0.
        zi = 0.
        delta=0.025
        rrx = random.normal(0., delta, natom)
        rry = random.normal(0., delta, natom)
        rrz = random.normal(0., delta, natom)
        with open("fcc.txt", "w") as f:
            for nx in arange(self.l1) :
                for ny in arange(self.l2) :
                    for nz in arange(self.l3) :
                        self.x[j] = xi + ax*nx + rrx[j]
                        self.y[j] = yi + ay*ny + rry[j]
                        self.z[j] = zi + az*nz + rrz[j]
                        f.write( "  %d   %8.3f   %8.3f   %8.3f \n" % (j, self.x[j], self.y[j], self.z[j]) )
                        j +=1
                        self.x[j] = xi + ax*nx + rrx[j] + 0.5*ax
                        self.y[j] = yi + ay*ny + rry[j] + 0.5*ay
                        self.z[j] = zi + az*nz + rrz[j]
                        f.write( "  %d   %8.3f   %8.3f   %8.3f \n" % (j, self.x[j], self.y[j], self.z[j]) )
                        j +=1
                        self.x[j] = xi + ax*nx + rrx[j] + 0.5*ax
                        self.y[j] = yi + ay*ny + rry[j]
                        self.z[j] = zi + az*nz + rrz[j] + 0.5*az
                        f.write( "  %d   %8.3f   %8.3f   %8.3f \n" % (j, self.x[j], self.y[j], self.z[j]) )
                        j +=1
                        self.x[j] = xi + ax*nx + rrx[j]
                        self.y[j] = yi + ay*ny + rry[j] + 0.5*ay
                        self.z[j] = zi + az*nz + rrz[j] + 0.5*az
                        f.write( "  %d   %8.3f   %8.3f    %8.3f \n" % (j, self.x[j], self.y[j], self.z[j]) )
                        j +=1
        print( "# end of initial fcc lattice construction")
    # Costruire file con il numero di vettori e le coordinate ridotte dei vettori
    # Modificare la routine fcc per leggere il file e costruire un cristallo date le coordinate dei vettori

    def gen_cr(self, mode, lx, ly, lz):
        from numpy import  random, arange, loadtxt, zeros
        # Unitary cell lengths
        ax = self.Lx/lx
        ay = self.Ly/ly
        az = self.Lz/lz
        print( "# VALUES BEFORE CELLS OPTIMIZATIONS")
        print( "# cell parameters a =%.4f  %.4f  %.4f" % (ax, ay, az) )
        print( "# mx = %d   my = %.d  mz = %.d" % (lx, ly, lz) )
        print( "# VALUES AFTER CELLS OPTIMIZATIONS")
        print( "# cell parameters a =%.4f  %.4f  %.4f" % (self.Lx/self.mx, self.Ly/self.my, self.Lz/self.mz))
        print( "# mx = %d   my = %.d  mz = %.d" % (self.mx,self.my, self.mz) )
        m = lx*ly*lz
        print( "# number of lattice cells m = %d" % m )
        j  = 0
        xi = 0.
        yi = 0.
        zi = 0.
        delta=0.001
        #   The file contains the vectors
        natom = self.nvec*m
        print( "# number of particles  %d " % natom )
        print( "# sides of md-box L = [ %.4f %.4f %.4f ]" % (self.Lx, self.Ly, self.Lz) )
        rrx = random.normal(0., delta, natom)
        rry = random.normal(0., delta, natom)
        rrz = random.normal(0., delta, natom)
        # Debugging (perfect cristal)
        # rrx = zeros(natom)
        # rry = zeros(natom)
        # rrz = zeros(natom)
        if mode > 0 :
        # at startup we randomly assign particle type
            self.Na = 0
            self.Nb = 0
            # prob = [self.na/N, 1.-self.na/N]
            # self.sp = random.choice(2,N,p=prob)
        with open("./files/input/cryst_in.txt", "w") as f:
            for nx in arange(lx) :
                for ny in arange(ly) :
                    for nz in arange(lz) :
                        for vec in arange(self.nvec) :
                            self.x[j] = xi + ax*nx + rrx[j] + self.xvector[vec]
                            self.y[j] = yi + ay*ny + rry[j] + self.yvector[vec]
                            self.z[j] = zi + az*nz + rrz[j] + self.zvector[vec]
                            self.sp[j]= self.spv[vec]
                            # print(self.sp[j])
                            f.write( "  %d   %8.3f   %8.3f   %8.3f   %d\n" % (j, self.x[j], self.y[j], self.z[j], self.sp[j]) )
                            if self.sp[j] == 0:
                                self.m[j] = self.ma
                                self.Na += 1
                            else:
                                self.m[j] = self.mb
                                self.Nb += 1
                            j +=1

        print( "# end of initial lattice construction")


    def write_gdr(self, N1, N2, T, gdr_out='gdr.out'):
    # """ T = number of times calcgdr was called """
        from numpy import zeros, pi, savetxt, column_stack
        V = zeros(self.kg)
        r = zeros(self.kg)
        g = zeros( (self.kg,3) )
        vol = self.Lx*self.Ly*self.Lz
        for lm in range(self.kg):
            r[lm] = (lm+0.5)*self.ldel
            V[lm] = 4./3.*pi*(self.ldel**3)*(3*lm*lm +3*lm + 1)
            g[lm,0] = self.gcount[lm,0]/(V[lm]*N1*(N1-1)*T/vol)
            g[lm,1] = self.gcount[lm,1]/(V[lm]*N1*N2*2*T/vol)
            g[lm,2] = self.gcount[lm,2]/(V[lm]*N2*(N2-1)*T/vol)
        gout = column_stack( (r, g[:,0], g[:,1], g[:,2]) )
        savetxt(gdr_out, gout , fmt=('%10.5g ','%12.7g','%12.7g ','%12.7g'), header="    'r'     'g11(r)'     'g12(r)'     'g22(r)'" )
        self.integrate_gdr(g)

    def write_integral(self, int_out='int.out'):
        # """ T = number of times calcgdr was called """
        from numpy import zeros, pi, savetxt, column_stack
        V = zeros(self.kg)
        r = zeros(self.kg)
        g = zeros( (self.kg,4) )
        # vol = self.Lx*self.Ly*self.Lz
        for lm in range(self.kg):
            r[lm] = (lm+0.5)*self.ldel
            V[lm] = 4./3.*pi*(self.ldel**3)*(3*lm*lm +3*lm + 1)
            g[lm,0] = self.integral[lm,0]
            g[lm,1] = self.integral[lm,1]
            g[lm,2] = self.integral[lm,2]
            g[lm,3] = self.integral[lm,3]
        intout = column_stack( (r, g[:,0], g[:,1], g[:,2], g[:,3]) )
        savetxt(int_out, intout , fmt=('%10.5g ','%12.7g','%12.7g ','%12.7g','%12.7g'), header="    'r'     'cord11(r)'     'cord12(r)'     'cord22(r)'     cord21(r)" )
