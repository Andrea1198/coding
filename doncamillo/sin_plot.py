import matplotlib.pyplot as plt
import numpy as np

file = "./files/d4opr_sin.dat"
data = np.loadtxt(file, unpack=True)
x = data[0,:]
ys = data[1:,:]
for y in ys:
    plt.plot(x, y)
plt.show()
