import numpy as np
import matplotlib.pyplot as plt
import random
import math as m
import csv

data = np.genfromtxt('realData.csv', delimiter=',')
print(data.shape)

step = data[:,0]
time = data[:,1]
position = data[:,2]
counter_t = data[:,3]
counter_c = data[:,4]


