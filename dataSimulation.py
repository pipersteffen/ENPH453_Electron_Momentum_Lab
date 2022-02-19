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

# For fermi distribution
HBAR = 6.582119569 * 10**(-16) #reduced plank's constant [eV⋅s]
E_MASS = 0.51099895000*1000000 #electron mass [MeV/c^2] -> [eV/c^2]

def electronType() -> str:
    return random.choice([0,1]) # 0=core, 1=valence

class fermiPDF(st.rv_continuous):
    def _pdf(self,E):
         densityOfStates = (1/(2*m.pi**2))*((2*E_MASS)/(HBAR**2)**(3/2))*E**(0.5)
         return densityOfStates

def energyToMomentum(e):
    return m.sqrt(e*2*E_MASS)

def initFermiEnergy(nFermiElements = 15000):
    #valid fermi energy range
    fermiRange = [0, energyToMomentum(7)]     # [ev]
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

def getRandFermiMomentum(energy, normProb):
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
            # assign momenta based on gaussian distribution, where sigma is converted to momentum from energy
            z_momenta[i] = random.gauss(0, energyToMomentum(4)) # mu=0, sigma=4 eV (energy)

        elif e_type[i] == 1: #'valence':
            total_momenta = getRandFermiMomentum(fermi_energy, fermi_energy_PDF)
            
            phi_random = random.randrange(0,180)
            # assign energy based on fermi momentum distribution
            z_momenta[i] = m.cos(m.radians(phi_random))*total_momenta

        else:
            print('Something went very, very wrong')

        # math in the OneNote, Tuesday Feb 8th notes
        theta[i] = m.atan(z_momenta[i] / PHOTON_ENERGY ) # 511 keV, guassian in eV

        # find height that particle will hit at
        heights[i] = m.tan(theta[i]) * DETECTOR_DIST

    ########################################################
    # Plots
    ########################################################
    
    # Histogram height distribution for each electron type
    core_heights = np.where(e_type == 0, heights, np.nan)
    valence_heights = np.where(e_type == 1, heights, np.nan)
    fig, (ax1, ax2) = plt.subplots(1,2)

    ax1.hist(core_heights, bins=numbins_height, alpha=0.8, color='skyblue', label='Core electrons', zorder=3)
    ax1.hist(valence_heights, bins=numbins_height, alpha=0.8, color='mediumseagreen', label='Valence electrons', zorder=3)
    ax1.legend()
    ax1.set_ylabel('Counts')
    ax1.set_xlabel('Height at detector plane [cm]')
    ax1.set_title('Impact height distribution at detector plane by electron type')
    ax1.grid(linestyle=':', zorder=0)

    # Histogram height distribution for each electron type
    core_mom = np.where(e_type == 0, z_momenta, np.nan)
    valence_mom = np.where(e_type == 1, z_momenta, np.nan)
    ax2.hist(core_mom/1000, bins=numbins_height, alpha=0.8, color='skyblue', label='Core electrons', zorder=3)
    ax2.hist(valence_mom/1000, bins=numbins_height, alpha=0.8, color='mediumseagreen', label='Valence electrons', zorder=3)
    ax2.legend()
    ax2.set_ylabel('Counts')
    ax2.set_xlabel('Momenta [keV]')
    ax2.set_title('Momentum distribution at detector plane by electron type')
    ax2.grid(linestyle=':', zorder=0)
    
    plt.show()


if __name__=="__main__":
    main()
