from dataclasses import replace
import numpy as np
def read(filename):
    data = []
    with open(filename, "r") as f:
        lines = f.readlines()
        for l in lines:
            l = l.split()
            while " " in l:
                l.remove(" ") 

            if l[0] != "#":
                for i in range(len(l)):
                    l[i] = float(l[i])
                data.append(l)
        
    return np.array(data)