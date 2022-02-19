from ast import Global
from cmath import pi
import numpy as np 
import random 
import matplotlib.pyplot as plt
import scipy.stats as st
import math as m 


class my_pdf(st.rv_continuous):
    def _pdf(self,E):
         #reduced plank's constant = 6.582119569...×10−16 [eV⋅s]
         hbar = 6.582119569 * 10**(-16)
         #electron mass = 0.51099895000 [MeV/c^2] 
         mass = 0.51099895000*1000000 # [eV]
         densityOfStates = (1/(2*m.pi**2))*((2*mass)/(hbar**2)**(3/2))*E**(0.5)
         return densityOfStates

mass = 0.51099895000*1000000 # [eV]
upper = (2*mass*7.0)**(0.5)

print(upper)
print((2*mass*4)**(0.5))
my_cv = my_pdf(a=0, b=upper, name='my_pdf')

numElements = 15000

#set energy range 
energy = np.linspace(0, upper, numElements)
energyKeV = energy/1000
#calculate DOS 
dos = my_cv.pdf(energy)

normalizer = 1/float(sum(dos))
normProb = [element * normalizer for element in dos]

#np.random.choice needs to have probability array sum to one 
#normalizing process is close to sum 1, but not exact as it is discrete, not continuous
#find difference from 1 and add to last element (as should have the highest probability)
diff = 1 - sum(normProb) 
normProb[len(normProb)-1] = normProb[len(normProb)-1] + diff

#test to ensure the same distribution is chosen 
N = 10000 
dist = np.empty(N) #random eV taken from Fermi distributino 

for i in range (0, N):
    dist[i] = np.random.choice(energy, size=1, p=normProb)

# show distribution for random choosing from Fermi momentum 
# overlay with the expected DOS ? 

# plt.figure()
# plt.plot(energyKeV, dos, color = 'black')
# plt.ylabel("DOS")
# plt.xlabel('Momentum [eV]')
# plt.show()

plt.figure()
plt.plot(energyKeV, normProb, color = 'black')
plt.ylabel("Normalized DOS")
plt.xlabel('Momentum [keV]')
plt.grid(linestyle = ':')
plt.title("Normalized Density of States vs. Momentum")
plt.show()

# plt.figure()
# plt.hist(dist, alpha = 0.7, bins = 80)
# plt.xlabel('Energy [eV]')
# plt.ylabel('Number of Occurrences')
# plt.show()
