import numpy as np
import matplotlib.pyplot as plt
import random
import math as m

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

def electronType() -> str:
    return random.choice([0,1]) # 0=core, 1=valence

def main():

    e_type = np.empty(N).astype('int')        # electron type
    z_energy = np.empty(N)      # energy in z direction
    theta = np.empty(N)         # angle
    heights = np.empty(N)        # height that particle will hit the detector

    for i in range(0, N):
        # Get random choice of valence or core electron hit by positron (step 3)
        e_type[i] = electronType()

        # Choose the z component of the momentum (step 4)
        if e_type[i] == 0: #'core':
            # assign energy based on gaussian distribution
            z_energy[i] = random.gauss(1, 4) # mu=1, sigma=4 eV

        elif e_type[i] == 1: #'valence':
            z_energy[i] = random.uniform(-10,10)
            # assign energy based on fermi momentum distribution
            # TODO

        else:
            print('Something went very, very wrong')

        # Calculate the angle??? (step 7 in pseudo-pseudo-code in OneNote)
        # math in the OneNote, Tuesday Feb 8th notes
        theta[i] = m.atan(z_energy[i] / PHOTON_ENERGY ) # 511 keV, guassian in eV

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


    plt.show()

    

if __name__=="__main__":
    main()