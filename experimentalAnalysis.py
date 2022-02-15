import numpy as np
import matplotlib.pyplot as plt
import math as m
from scipy.optimize import curve_fit
from scipy.stats import norm

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

    #counts = counter_coinc/counter_total # normalize hits on detector
    #heights = counter_coinc**2/counter_total

    counts = counter_coinc**2/counter_total # corrected # of counts
    err_counts = ((2*counter_coinc/counter_total * counter_coinc**0.5)**2 + (counter_coinc**2/counter_total**2 * counter_total**0.5)**2)**0.5 # uncertainty in counts

    x = position*0.5*STEPSIZE
    n = len(x)
    mean = sum(x*counts)/n
    sigma = m.sqrt(sum(counts*(x-mean)**2)/n)

    # fit data to generic gaussian function
    popt, pcov = curve_fit(gauss, x, counts, p0=[max(counts),mean,sigma], sigma=err_counts) # popt [a, x0, sigma]
    
    offset = popt[1]
    #print('Curve fit x0 = ', x0) # -0.06188926436973738

    if FIT_TAILS:
        # Fit tails of data to gaussian function
        # Need subset of data: furthest left  - left+ 0.793, furthest right-0.793 - furthest right
        #delta = 0.793
        delta = 0.45 #arbitrary, eyeballed from data

        dL = offset - delta # left side
        dR = offset + delta
        
        # subselect the x values of points on the tails
        left_points = x[x <= dL]
        right_points = x[x >= dR]
        tail_x = np.concatenate([left_points, right_points])

        mid_x = x[(x >= dL) & (x <= dR)]

        # subselect the count values of points on the tails
        tail_counts_left = counts[x <= dL]
        tail_counts_right = counts[x >= dR]
        tail_counts = np.concatenate([tail_counts_left, tail_counts_right])

        mid_counts = counts[(x >= dL) & (x <= dR)]
        #mid_counts = np.where(np.logical_and(x >= -delta, x <= delta))
        print(mid_counts)
        #mid_counts = np.concatenate([counts[x >= -delta], counts[x <= delta]])

        # subselect the error counts array to points on the tails
        err_counts_left = err_counts[x <= dL]
        err_counts_right = err_counts[x >= dR]
        err_counts_tail = np.concatenate([err_counts_left, err_counts_right])

        popt2, pcov2 = curve_fit(gauss, tail_x, tail_counts, p0=[max(counts),mean,0.3], sigma=err_counts_tail) # 0.3 is a guess from the data points dist

        tail_fit = gauss(x, *popt2)
        all_fit = gauss(x, *popt)

        # plot the fit
        fig3, ax3 = plt.subplots()
        ax3.plot(x, counts,  '.k', label='All data')
        ax3.plot(x, all_fit, ':k', label='All data fit')
        ax3.plot(tail_x, tail_counts, '.r', label='Tail data')
        ax3.plot(x, tail_fit, '-r', label='Tail data fit')

        # Now fit the parabola part....
        # # subtract tail fit values from data TODO
        # tail_fit_mid = np.concatenate([tail_fit[x >= -delta], tail_fit[x <= delta]])
        mid_tail_fit = tail_fit[(x >= dL) & (x <= dR)]
        counts_sub = mid_counts - mid_tail_fit

        # fit quadtratic to that
        popt3, pcov3 = curve_fit(poly, mid_x, counts_sub)#, p0=[-8, -1, 1.5]) # need a,b,c guess for p0

        poly_fit = poly(mid_x, *popt3)
        ax3.plot(mid_x, poly_fit, '-b', label='mid subtract fit')
        ax3.plot(mid_x, counts_sub, '*', label='subtracted')

        intercept = np.roots(popt3)
        print(intercept)
        print('fermi vals: ')
        print(intercept[0]-offset)
        print(intercept[1]-offset)

        ax3.set_xlabel('Position on detector plane')
        ax3.set_ylabel('Corrected number of hits on detector')
        ax3.set_title('Experimental Data')
        ax3.legend()

        # FIND RATIO
        poly_sum = np.sum(poly_fit[poly_fit>=0])
        tail_sum = np.sum(tail_fit)
        all_sum = np.sum(all_fit)

        print('poly coeffs:',popt3)

        # percent of count is tail/core
        perc_core = tail_sum/(tail_sum+poly_sum)

        print(poly_fit)

        # percent of count is poly/valence
        perc_val = poly_sum/(tail_sum+poly_sum)
        
        print('Percent core: ', perc_core*100)
        print('Percent valence: ', perc_val*100)

        print('ratio or poly curve to tail: ', poly_sum/tail_sum)


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
        ax2.bar(x, counts, width=0.02)
        ax2.set_xlabel('Position on detector plane')
        ax2.set_ylabel('Corrected counts per total coincidence hits')
        ax2.set_title('Experimental Data')
    
    if PLOT_DATA_FIT or PLOT_DATA_BAR or PLOT_TAIL_FIT:
        plt.show()

if __name__=="__main__": 
    main()