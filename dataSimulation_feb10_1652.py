import numpy as np
import matplotlib.pyplot as plt
import random
import math as m
import scipy.stats as st

# Number of times to run
N = 10000
STEPSIZE = 0.0635   # cm
HMAX = 2.667        # cm, 1 degree above or below axis

# num of elements: 85
# num = m.ceil(2*HMAX/STEPSIZE + 1)
# print(num)
# exit()
numbins_height = 85
HEIGHTS = np.linspace(-HMAX, HMAX, numbins_height)  # use for binning, needs to go through 0

PHOTON_ENERGY = 511000 # eV
DETECTOR_DIST = 150 # cm

#for fermi distribution
HBAR = 6.582119569 * 10**(-16) #reduced plank's constant [eV⋅s]
E_MASS = 0.51099895000*1000000 #electron mass [MeV/c^2] -> [eV/c^2]

def electronType() -> str:
    return random.choice([0,1]) # 0=core, 1=valence

class fermiPDF(st.rv_continuous):
    def _pdf(self,E):
         densityOfStates = (1/(2*m.pi**2))*((2*E_MASS)/(HBAR**2)**(3/2))*E**(0.5)
         return densityOfStates

def initFermiEnergy(nFermiElements = 15000):
    #valid fermi energy range
    fermiRange = [0, 7]     # [ev]
    # valid fermi energies [eV]
    fermiEnergy = np.linspace(fermiRange[0], fermiRange[1], nFermiElements)

    #initialize fermi energy distribution and density of states 
    fermi = fermiPDF(a = fermiRange[0], b = fermiRange[1], name='fermiPDF')
    dos = fermi.pdf(fermiEnergy)

    #normalize density of states so integral = 1
    normalizer = 1/float(sum(dos))
    normProb = [element * normalizer for element in dos]

    #small correction factor for discrete integral 
    diff = 1 - sum(normProb) 
    normProb[len(normProb)-1] = normProb[len(normProb)-1] + diff

    return fermiEnergy, normProb

def getRandFermiEnergy(energy, normProb):
    return np.random.choice(energy, size=1, p=normProb)


def main():

    e_type = np.empty(N).astype('int')        # electron type
    z_energy = np.empty(N)                    # energy in z direction
    z_momenta = np.empty(N)                   # momenta in z direction 
    theta = np.empty(N)                       # angle
    heights = np.empty(N)                     # height that particle will hit the detector

    #initialize fermi energy distribution
    fermi_energy, fermi_energy_PDF = initFermiEnergy()

    for i in range(0, N):
        # Get random choice of valence or core electron hit by positron (step 3)
        e_type[i] = electronType()

        # Choose the z component of the momentum (step 4)
        if e_type[i] == 0: #'core':
            # assign energy based on gaussian distribution
            z_energy[i] = random.gauss(1, 4) # mu=1, sigma=4 eV
            if z_energy[i] < 0:
                z_momenta[i] = (-1)*(abs(z_energy[i])*2*E_MASS)**(1/2)
            else: 
                z_momenta[i] = (z_energy[i]*2*E_MASS)**(1/2)

        elif e_type[i] == 1: #'valence':
            total_fermi = getRandFermiEnergy(fermi_energy, fermi_energy_PDF)
            total_momenta = (total_fermi*2*E_MASS)**(1/2)
            
            phi_random = random.randrange(0,180)
            # assign energy based on fermi momentum distribution
            z_energy[i] = m.cos(m.radians(phi_random))*total_fermi
            z_momenta[i] = m.cos(m.radians(phi_random))*total_momenta

        else:
            print('Something went very, very wrong')

        # math in the OneNote, Tuesday Feb 8th notes
        theta[i] = m.atan(z_momenta[i] / PHOTON_ENERGY ) # 511 keV, guassian in eV

        # find height that particle will hit at
        heights[i] = m.tan(theta[i]) * DETECTOR_DIST

    ########################################################
    # Plot stuff
    ########################################################

    # Histogram the z energy distribution for each electron type
    core_energy = np.where(e_type == 0, z_energy, np.nan)
    valence_energy = np.where(e_type == 1, z_energy, np.nan)
    plt.figure(1)
    plt.hist(core_energy, bins=numbins_height, alpha=0.5, label='core')
    plt.hist(valence_energy, bins=numbins_height, alpha=0.5, label='valence')
    plt.xlabel('Z energy [eV]')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('Distribution of z-component of energies by electron type')
    
    # Histogram height distribution for each electron type
    core_heights = np.where(e_type == 0, heights, np.nan)
    valence_heights = np.where(e_type == 1, heights, np.nan)
    plt.figure(2)
    plt.hist(core_heights, bins=numbins_height, alpha=0.5, label='core')
    plt.hist(valence_heights, bins=numbins_height, alpha=0.5, label='valence')
    plt.legend()
    plt.ylabel('Counts')
    plt.xlabel('Height at detector plane [cm]')
    plt.title('Distribution of height of particle at detector plane by electron type')

    # Histogram height distribution for each electron type
    core_mom = np.where(e_type == 0, z_momenta, np.nan)
    valence_mom = np.where(e_type == 1, z_momenta, np.nan)
    plt.figure(3)
    plt.hist(core_mom, bins=numbins_height, alpha=0.5, label='core')
    plt.hist(valence_mom, bins=numbins_height, alpha=0.5, label='valence')
    plt.legend()
    plt.ylabel('Counts')
    plt.xlabel('Momenta [eV?]')
    plt.title('Distribution of momenta of particle at detector plane by electron type')
    plt.show()

    

if __name__=="__main__":
    main()