import numpy as np
import matplotlib.pyplot as plt
import math as m
from scipy.optimize import curve_fit

STEPSIZE = 0.0635 # cm

PLOT_DATA_FIT = 0
PLOT_DATA_BAR = 0
FIT_TAILS = 1
PLOT_TAIL_FIT = 1

def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def poly(x, a, b, c):
    # ax**2 + bx + c
    return a*x**2 + b*x + c

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

    if FIT_TAILS:
        # Fit tails of data to gaussian function
        # Need subset of data: furthest left  - left+ 0.793, furthest right-0.793 - furthest right
        #delta = 0.793
        delta = 0.6 #arbitrary, eyeballed from data
        
        # subselect the x values of points on the tails
        left_points = x[x <= -delta]
        right_points = x[x >= delta]
        tail_x = np.concatenate([left_points, right_points])

        # subselect the count values of points on the tails
        tail_counts_left = counts[x <= -delta]
        tail_counts_right = counts[x >= delta]
        tail_counts = np.concatenate([tail_counts_left, tail_counts_right])
        
        print('tail x: ', tail_x)
        print('tail counts: ', tail_counts)

        # fit the gauss
        n2 = len(tail_x)
        mean2 = sum(tail_x*tail_counts)/n2
        sigma2 = m.sqrt(sum(tail_counts*(tail_x - mean2)**2)/n2)

        popt2, pcov2 = curve_fit(gauss, tail_x, tail_counts, p0=[max(tail_counts),mean2,sigma2]) # popt [a, x0, sigma]

        print('popt2:', popt2)
        print('pcov2:', pcov2)
    
        # plot the fit
        fig3, ax3 = plt.subplots()
        ax3.plot(x, counts,  '.', label='all data')
        ax3.plot(tail_x, tail_counts, '.', label='tail points')
        ax3.plot(tail_x, gauss(tail_x, *popt2), label='tail fit')
        ax3.set_xlabel('Position on detector plane')
        ax3.set_ylabel('Number of hits on detector')
        ax3.set_title('Experimental Data')
        ax3.legend()

        # Now fit the parabola part....
  
        # subtract tail fit values from data TODO
        # x_prime = 

        # fit parabola TODO
        # popt, pcov = curve_fit(poly, x, counts)

        # plot the fit
       

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
    
    if PLOT_DATA_FIT or PLOT_DATA_BAR or PLOT_TAIL_FIT:
        plt.show()

if __name__=="__main__": 
    main()