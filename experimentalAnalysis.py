import numpy as np
import matplotlib.pyplot as plt
import math as m
import uncertainties as unc
from uncertainties import unumpy
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

def quadraticEquation(fit):
    a = fit[0]
    b = fit[1]
    c = fit[2]
    d = b**2 - 4*a*c
    sol1 = (-b - d**0.5)/(2*a)
    sol2 = (-b + d**0.5)/(2*a)
    return [sol1, sol2]

def main(): 
    data = np.genfromtxt('realData.csv', delimiter=',', skip_header=True)
    position = data[:,2]
    counter_total = data[:,3]
    counter_coinc = data[:,4]

    counts = counter_coinc**2/counter_total # corrected # of counts
    err_counts = ((2*counter_coinc/counter_total * counter_coinc**0.5)**2 + (counter_coinc**2/counter_total**2 * counter_total**0.5)**2)**0.5 # uncertainty in counts

    x = position*0.5*STEPSIZE
    n = len(x)
    mean = sum(x*counts)/n
    sigma = m.sqrt(sum(counts*(x-mean)**2)/n)

    # fit data to generic gaussian function
    popt, pcov = curve_fit(gauss, x, counts, p0=[max(counts),mean,sigma], sigma=err_counts) # popt [a, x0, sigma]
    
    offset = unc.ufloat(popt[1], np.sqrt(pcov[1,1]))
    #print('Curve fit x0 = ', x0) # -0.06188926436973738

    if FIT_TAILS:
        # Fit tails of data to gaussian function
        # Need subset of data: furthest left  - left+ 0.793, furthest right-0.793 - furthest right
        #delta = 0.793
        delta = 0.45 #arbitrary, eyeballed from data

        dL = offset - delta # left side
        dR = offset + delta
        
        # Subselect the x values of points on the tails
        left_points = x[x <= dL]
        right_points = x[x >= dR]
        tail_x = np.concatenate([left_points, right_points])

        mid_x = x[(x >= dL) & (x <= dR)]

        # Subselect the count values of points on the tails
        tail_counts_left = counts[x <= dL]
        tail_counts_right = counts[x >= dR]
        tail_counts = np.concatenate([tail_counts_left, tail_counts_right])

        mid_counts = counts[(x >= dL) & (x <= dR)]

        # Subselect the error counts array to points on the tails
        err_counts_left = err_counts[x <= dL]
        err_counts_right = err_counts[x >= dR]
        err_counts_tail = np.concatenate([err_counts_left, err_counts_right])

        popt2, pcov2 = curve_fit(gauss, tail_x, tail_counts, p0=[max(counts),mean,0.3], sigma=err_counts_tail) # 0.3 is a guess from the data points dist

        tail_fit = gauss(x, *popt2)
        all_fit = gauss(x, *popt)

        # Now fit the parabola part
        mid_tail_fit = tail_fit[(x >= dL) & (x <= dR)]
        counts_sub = mid_counts - mid_tail_fit
        err_counts_sub = abs(counts_sub)**0.5

        # fit quadratic to that, with errors
        popt3Temp, pcov3Temp = curve_fit(poly, mid_x, counts_sub)
        popt3, pcov3 = curve_fit(poly, mid_x, counts_sub, sigma=err_counts[(x >= dL) & (x <= dR)], p0=popt3Temp)
        quadFitError = unumpy.uarray(popt3, np.sqrt(np.diag(pcov3)))

        poly_fit = poly(mid_x, *popt3)

        intercept = quadraticEquation(quadFitError)
        average_int = (abs(intercept[0] - offset) + abs(intercept[1] - offset))/2
        fermiMomentum = average_int*511/150 #Momentum in keV

        # FIND RATIO OF CORE TO VALENCE ELECTRONS
        poly_sum = np.sum(poly_fit[poly_fit>=0])
        poly_err = (poly_sum)**0.5
        tail_sum = np.sum(tail_fit)
        tail_err = tail_sum**0.5
        polyUnc = unc.ufloat(poly_sum, poly_err)
        tailUnc = unc.ufloat(tail_sum, tail_err)

        # percent of count is tail/core
        perc_core = tailUnc/(tailUnc+polyUnc)
        # percent of count is poly/valence
        perc_val = polyUnc/(tailUnc+polyUnc)
        
        # generating errors in each fit coefficients
        allFitError = unumpy.uarray(popt, np.sqrt(np.diag(pcov)))
        tailFitError = unumpy.uarray(popt2, np.sqrt(np.diag(pcov2)))

        print('Percent core: ', perc_core*100)
        print('Percent valence: ', perc_val*100)
        print('Ratio or poly curve to tail: ', poly_sum/tail_sum)
        print('Fermi Momentum of fitted distrubtion: ', fermiMomentum, " keV.")
        print('Error in All Fit: ', allFitError)
        print('Error in Tail Fit: ', tailFitError)
        print('Error in Quadratic Fit: ', quadFitError)


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

    if PLOT_TAIL_FIT:
        fig3 = plt.figure(figsize=(16,10))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(2, 2, 2)
        ax3 = plt.subplot(2, 2, 4)

        fig3.suptitle('Experimental Data with Curve Fitting')
        ax1.errorbar(x, counts, yerr=err_counts, fmt='.k', capsize=3,label='Experimental Data')
        ax1.plot(x, all_fit, '--k', label='Gaussian Fit to All Data')
        ax1.set_xlabel('Position on Detector Plane')
        ax1.set_ylabel('Corrected number of hits on detector')
        ax1.legend()
        ax1.set_title("Gaussian Fit to All Data")
        ax1.grid(linestyle=":")

        ax2.errorbar(x, counts, yerr=err_counts, fmt='.k', capsize=3,label='Experimental Data')
        ax2.errorbar(tail_x, tail_counts, yerr=err_counts_tail, fmt='*r', capsize=3, label='Tail Data for Gaussian Fit')
        ax2.plot(x, tail_fit, '-r', label='Gaussian Fit to Tail Data')
        ax2.set_xlabel('Position on Detector Plane')
        ax2.set_ylabel('Corrected number of hits on detector')
        ax2.legend()
        ax2.set_title("Gaussian Fit to Tail Data")
        ax2.grid(linestyle=":")

        ax3.errorbar(x, counts, yerr=err_counts, fmt='.k', capsize=3,label='Experimental Data')
        ax3.errorbar(mid_x, counts_sub, yerr=err_counts[(x >= dL) & (x <= dR)], fmt='*b', capsize=3, label='Experimental Data Minus Tail Fit')
        ax3.plot(mid_x, poly_fit, '-b', label='Quadratic Fit to Subtracted Data')
        ax3.set_xlabel('Position on Detector Plane')
        ax3.set_ylabel('Corrected number of hits on detector')
        ax3.legend()
        ax3.set_title("Quadratic Fit to Centre Data")
        ax3.grid(linestyle=":")


    
    if PLOT_DATA_FIT or PLOT_DATA_BAR or PLOT_TAIL_FIT:
        plt.show()


if __name__=="__main__": 
    main()