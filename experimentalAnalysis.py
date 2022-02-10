import numpy as np
import matplotlib.pyplot as plt
import random
import math as m
import csv

STEPSIZE = 0.0635 # cm

data = np.genfromtxt('realData.csv', delimiter=',')
print(data.shape)

step = data[:,0]
time = data[:,1]
position = data[:,2]
counter_total = data[:,3]
counter_coinc = data[:,4]

counts = counter_coinc/counter_total

# Plot
# y value: counter_c / counter_total --> gaussian fit
# x value: postition * 0.5 * step size
    # in sim we get a coincidence hit every time we send a positron

plt.figure()
plt.hist(counter_coinc**2/counter_total, bins=51)
plt.show()