#Librerie
import math
import pandas
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def main():
    #Esportazione dei dati dai file e creazione degli array
    dati = pandas.read_csv("./gaas.dat")
    length = np.array(dati['Length(nm)'])
    n = np.array(dati['n'])

    #costanti
    h   = 1.05457182*1e-34
    m_e = 9.1093837*1e-31
    pi  = np.pi

    E_l = ( (h**2) * (pi**2) ) / 0.134 / m_e / length**2
    E_n = ( (h**2) * (pi**2) * n**2 ) / 0.134e-8 / m_e


    print(E_l)
    #grafico parte reale
    plt.plot(length, E_l, color="r", linewidth=2.5)
    plt.title("E vs L")
    plt.xlabel("L")
    plt.ylabel("E")
    plt.savefig("./E_l.png")
    plt.close()

    print(E_n)
    #grafico parte reale
    plt.plot(length, E_n, color="r", linewidth=2.5)
    plt.title("E vs n")
    plt.xlabel("n")
    plt.ylabel("E")
    plt.savefig("./E_n.png")
    plt.close()

if __name__ == "__main__":
    main()
