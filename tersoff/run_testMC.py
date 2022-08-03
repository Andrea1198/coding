"""
Units of the system:
    Mass: uma (1.661e-27 Kg)
    Time: 1.0364e-14 s (~10 fs)
    Length: Angstrom (1e-10 m)
    Temperature: 1e5 K
    
How to convert temperature to energy:
kb in this units is (J/K -> Kg*m/s/s/K)
"""


from numpy import zeros,pi,sin,cos,empty,array,float64,random,int32
from sympy import Eq
from src.ter import TER
from time import process_time
from src.eqmc import eqmc
from src.read_steps import read
if __name__ == "__main__":
    # use a fixed seed for debugging purposes
    # random.seed(1234567)
    dir = "./files/output/"
    # r-parameters # defaults
    count= 0
    freq = 1       
    mx   = 4       
    my   = 3       
    mz   = 4       
    dt   = 0.05    
    mode = 2        # random particle exchanges & initial velocities
    Na   = 0        # initial fraction of type a particles (< N)
    delta= 0.01     # Random position movement normalized to 1
    dV0  = 0.01      # Random volume change normalized to 1
    p_ext= .097
    
    #choice = (input("Insert 0 to run 'steps.txt' as input file, or 1 to write steps yourself\n"))
    choice = 0
    #if choice != "0" or choice != "1": choice = 0
    #else: choice = int(choice)
    if choice == 0:
        nstep, kt, freq, nstep_gdr, freq_gdr = read("./files/input/steps.txt").transpose()
        nstep   = int32(nstep)
        freq    = int32(freq)
        nstep_gdr    = int32(nstep_gdr)
        freq_gdr    = int32(freq_gdr)
        n_pas = len(nstep)
        filename = [""] * n_pas
        filename[0]     = "crystal"
        filename[2]     = "liquid"
        filename[-1]    = "glass"

        print( "# integration time step dt = %8.4f" % dt )
        EQ=TER(mx, my, mz)
        N=EQ.npart
        # rhofact= EQ.lb/(EQ.Lx*EQ.Ly*EQ.Lz)
        print( "# sides of the rectangular MD-box L = [ %8.4f %8.4f %8.4f ]" % (EQ.Lx, EQ.Ly, EQ.Lz) )
        print( "# density rho = %8.4f"  % (EQ.rho ) )
        tic = process_time()
        EQ.gen_cr(1, mx, my, mz)
        print( "# number of particles of %s  N(%s) = %d" % ('O ', 'O ', EQ.Na) )
        print( "# number of particles of %s  N(%s) = %d" % ('Si', 'Si', EQ.Nb) )
        EQ.writexyz(N, dir+'fcc.xyz')
        EQ.write_input(N, 0, conf_out=dir+'conf_in.csv')
        print("# cpu time for cristal startup tau=%12.5f" % (process_time()-tic))
        tic = process_time()
        estep=EQ.read_input(N, conf_in=dir+'conf_in.csv')
        print("# cpu time for reading initial configuration tau=%12.5f" % (process_time()-tic))
        tic = process_time()

        for i in range(n_pas):
            if nstep_gdr[i] != 0: temp = ""
            else: temp = "o"
            print("# Step {}, kt = {} °K w/{} gdr".format(i, kt[i]*1.e4, temp))
            # print( "# mean (kinetic) temperature kt = %8.4f °K" % (kt[i]*1.e4) )
            eqmc(EQ, N, kt[i], nstep[i], dt, freq[i], 2, dir, delta, dV0, p_ext)
            if nstep_gdr[i] != 0:
                # print( "# Computing gdr ... " )
                eqmc(EQ, N, kt[i], nstep_gdr[i], dt, freq_gdr[i], 1, dir, delta, dV0, p_ext)
                if filename[i] != "":
                    EQ.write_gdr(EQ.Na, EQ.Nb, nstep_gdr[i]/freq_gdr[i], gdr_out=dir+'gdr_{}.out'.format(filename[i]))
                    EQ.write_integral(int_out=dir+ 'int_{}.out'.format(filename[i]))
            if filename[i] != "":
                EQ.write_input(N, 0, conf_out=dir+filename[i] + '.csv')
                EQ.writexyz(N, dir+"conf_{}.xyz".format(filename[i]))
            EQ.gcount = zeros((EQ.kg,3))
            EQ.integral = zeros((EQ.kg,4))
            count += nstep[i] + nstep_gdr[i]

            print("# Final length of the box: {} {} {}".format(EQ.Lx, EQ.Ly, EQ.Lz))
    else:
        
        kt, nstep, freq, nstep_gdr, freq_gdr = input("Insert: kt, nsteps, freq, nstep_gdr, freq_gdr (0 to stop)\n").split()
        # kt = int(kt)
        # nstep = int(nstep)
        # freq = int(freq)
        # nstep_gdr = int(nstep_gdr)
        # freq_gdr = int(freq_gdr)        
        if nstep_gdr*nstep != 0: 
            filename = input("Insert a filename\n")
        count = 0
        print( "# integration time step dt = %8.4f" % dt )
        EQ=TER(mx, my, mz)
        N=EQ.npart
        # rhofact= EQ.lb/(EQ.Lx*EQ.Ly*EQ.Lz)
        print( "# sides of the rectangular MD-box L = [ %8.4f %8.4f %8.4f ]" % (EQ.Lx, EQ.Ly, EQ.Lz) )
        print( "# density rho = %8.4f"  % (EQ.rho ) )
        tic = process_time()
        EQ.gen_cr(1, mx, my, mz)
        print( "# number of particles of %s  N(%s) = %d" % ('O ', 'O ', EQ.Na) )
        print( "# number of particles of %s  N(%s) = %d" % ('Si', 'Si', EQ.Nb) )
        EQ.writexyz(N, dir+'fcc.xyz')
        EQ.write_input(N, 0, conf_out=dir+'conf_in.csv')
        print("# cpu time for cristal startup tau=%12.5f" % (process_time()-tic))
        tic = process_time()
        estep=EQ.read_input(N, conf_in=dir+'conf_in.csv')
        print("# cpu time for reading initial configuration tau=%12.5f" % (process_time()-tic))
        tic = process_time()
        while nstep != 0:
            nstep       = int32(nstep)
            freq        = int32(freq)
            nstep_gdr   = int32(nstep_gdr)
            freq_gdr    = int32(freq_gdr)
            kt          = float64(kt)

            if nstep_gdr != 0: temp = ""
            else: temp = "o"
            print("# Step {}, kt = {} °K w/{} gdr".format(count, kt*1.e4, temp))

            # print( "# mean (kinetic) temperature kt = %8.4f °K" % (kt[i]*1.e4) )
            eqmc(EQ, N, kt, nstep, dt, freq, 2, dir, delta, dV0, p_ext)
            if nstep_gdr != 0:
                # print( "# Computing gdr ... " )
                eqmc(EQ, N, kt, nstep_gdr, dt, freq_gdr, 1, dir, delta, dV0, p_ext)
                if filename != "":
                    EQ.write_gdr(EQ.Na, EQ.Nb, nstep_gdr/freq_gdr, gdr_out=dir+'gdr_{}.out'.format(filename))
                    EQ.write_integral(int_out=dir+ 'int_{}.out'.format(filename))
            if filename != "":
                EQ.write_input(N, 0, conf_out=dir+filename + '.csv')
            EQ.gcount = zeros((EQ.kg,3))
            EQ.integral = zeros((EQ.kg,4))
            count += nstep + nstep_gdr

            print("# Final length of the box: {} {} {}".format(EQ.Lx, EQ.Ly, EQ.Lz))
            data_in = input("Insert: kt, nsteps, freq, nstep_gdr, freq_gdr (0 to stop)\n").split()
            if len(data_in) < 5: break
            kt = int(kt)
            nstep = int(nstep)
            freq = int(freq)
            nstep_gdr = int(nstep_gdr)
            freq_gdr = int(freq_gdr)
            if nstep_gdr * nstep != 0: 
                filename = input("Insert a filename\n")
            count += 1
