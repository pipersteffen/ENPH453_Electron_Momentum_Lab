import numpy as np
import matplotlib.pyplot as plt
import random
import math as m
import csv
from scipy.optimize import curve_fit
from scipy import exp
from sympy import N

STEPSIZE = 0.0635 # cm

def gauss(x, a, x0, sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def main(): 
    data = np.genfromtxt('realData.csv', delimiter=',', skip_header=True)
    print(data.shape)

    step = data[:,0]
    time = data[:,1]
    position = data[:,2]
    counter_total = data[:,3]
    counter_coinc = data[:,4]

    counts = counter_coinc/counter_total

    print(counts)
    x = position*0.5*STEPSIZE

    n = len(x)
    mean = sum(x*counts)/n
    sigma = m.sqrt(sum(counts*(x-mean)**2)/n)

    # Plot
    # y value: counter_c / counter_total --> gaussian fit
    popt, pcov = curve_fit(gauss, x, counts, p0=[max(counts),mean,sigma])
    # x value: postition * 0.5 * step size
        # in sim we get a coincidence hit every time we send a positron

    plt.figure(1)
    plt.plot(x, counts,  '.', label='data')
    plt.plot(x, gauss(x, *popt), label='fit')
    plt.xlabel('Position on detector plane')
    plt.ylabel('Number of hits on detector')
    plt.title('Experimental Data')
    plt.legend()
    plt.show()

if __name__=="__main__": 
    main()