import numpy as np
import matplotlib.pyplot as plt
import os


files_d4    = []
files_d1    = []
files_d4    = ["d4opr_gaus.dat", "d4opr_r4.dat", "d4opr_r5.dat", "d4opr_rm1.dat", "d4opr_sin.dat"] 
xlims_d4    = [[ 0., 20], [0,20], [0,20], [0,20], [0,20]]
ylims_d4    = [[-0.5, 2.5], [-3,3], [-3,3], [-1, 1.5], [-1.5, 1.5]]
labels_d4   = [["gaus", "r^4", "r^5", "r^-1", "sin"], ["4th der with spline"]*5]

for i,file in enumerate(files_d4):
    data = np.loadtxt("./files/"+file, unpack=True)
    r = data[0]
    ys= data[1:]
    plt.close()
    plt.xlim(xlims_d4[i])
    plt.ylim(ylims_d4[i])
    plt.scatter(r, ys[1], s=0.5, label=labels_d4[0][i])
    plt.scatter(r, ys[0], s=0.5, label=labels_d4[1][i])
    plt.legend()
    plt.savefig("./images/" + file[:-3] + "jpg")

files_d1    = ["dopr_gaus.dat", "dopr_r2.dat", "dopr_rm2.dat", "dopr_sin.dat"]
xlims_d1    = [[ 0., 20], [0,3], [0,20], [0,20]]
ylims_d1    = [[-1.5,1.5], [0,9], [-3,3], [-1.5,1.5]]
labels_d1   = [["gaus", "r^2", "r^-2", "sin"], ["1st der analitic"]*4, ["1st der bspline"]*4]

for i,file in enumerate(files_d1):
    data = np.loadtxt("./files/"+file, unpack=True)
    r = data[0]
    ys= data[1:]
    plt.close()
    plt.xlim(xlims_d1[i])
    plt.ylim(ylims_d1[i])
    plt.xlabel('r', fontsize=15)
    plt.ylabel("f(r), f'(r)", fontsize=15)

    plt.plot(r, ys[0], color="blue", linestyle="-.", label=labels_d1[0][i])
    plt.scatter(r, ys[1], s=10, color="black", label=labels_d1[1][i])
    plt.scatter(r, ys[2], s=1, color="red", label=labels_d1[2][i])
    plt.legend(markerscale=3, prop={'size' : 16})
    plt.savefig("./images/" + file[:-3] + "jpg")

data = np.loadtxt("./files/dopr_sin.dat", unpack=True)
r = data[0]
ys= data[1:]
plt.close()
plt.xlim(xlims_d1[3])
plt.ylim(ylims_d1[3])
plt.xlabel('r', fontsize=15)
plt.ylabel("f(r), f'(r)", fontsize=15)

plt.plot(r, ys[0], color="blue", linestyle="-.", label=labels_d1[0][3])
#plt.plot(r[::3], ys[1][::3], color="black", linestyle="-", label=labels_d1[1][i])
plt.scatter(r, ys[1], s=10, color="black", label=labels_d1[1][3])
#plt.plot(r, ys[2], color="red", linestyle="--", label=labels_d1[2][i])
plt.scatter(r, ys[2], s=1, color="red", label=labels_d1[2][3])
plt.legend(markerscale=3, prop={'size' : 16}, loc=2)
plt.savefig("./images/dopr_sin.jpg")

