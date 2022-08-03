import matplotlib.pyplot as plt
import numpy as np

file = "./files/dvpot_mat.dat"
data = np.loadtxt(file, unpack=True)
x = data[0,:]
ys = data[2:,:]
plt.ylim([-2,2])

for y in ys:
    plt.scatter(x, y, s=0.5)
plt.show()
