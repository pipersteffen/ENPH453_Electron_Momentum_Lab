import numpy as np
import matplotlib.pyplot as plt
import math as m
from scipy.optimize import curve_fit

STEPSIZE = 0.0635 # cm

PLOT_DATA_FIT = 0
PLOT_DATA_BAR = 0

def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def main(): 
    data = np.genfromtxt('realData.csv', delimiter=',', skip_header=True)

    step = data[:,0]
    time = data[:,1]
    position = data[:,2]
    counter_total = data[:,3]
    counter_coinc = data[:,4]

    counts = counter_coinc/counter_total # normalize hits on detector
    heights = counter_coinc**2/counter_total

    x = position*0.5*STEPSIZE
    n = len(x)
    mean = sum(x*counts)/n
    sigma = m.sqrt(sum(counts*(x-mean)**2)/n)

    # fit data to generic gaussian function
    popt, pcov = curve_fit(gauss, x, counts, p0=[max(counts),mean,sigma]) # popt [a, x0, sigma]
    
    x0 = popt[1]
    print('Curve fit x0 = ', x0) # -0.06188926436973738

    # Plot data & fit
    if PLOT_DATA_FIT:
        fig1, ax1 = plt.subplots()
        ax1.plot(x, counts,  '.', label='data')
        ax1.plot(x, gauss(x, *popt), label='fit')
        ax1.set_xlabel('Position on detector plane')
        ax1.set_ylabel('Number of hits on detector')
        ax1.set_title('Experimental Data')
        ax1.legend()

    if PLOT_DATA_BAR:
        fig2, ax2 = plt.subplots()
        ax2.bar(x, heights, width=0.02)
        ax2.set_xlabel('Position on detector plane')
        ax2.set_ylabel('Corrected counts per total coincidence hits')
        ax2.set_title('Experimental Data')
    
    if PLOT_DATA_FIT or PLOT_DATA_BAR:
        plt.show()

if __name__=="__main__": 
    main()