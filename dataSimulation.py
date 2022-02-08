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

    e_type = np.empty(N)        # electron type
    z_energy = np.empty(N)      # energy in z direction
    theta = np.empty(N)         # angle
    height = np.empty(N)        # height that particle will hit the detector

    for i in range(0, N):
        # Get random choice of valence or core electron hit by positron (step 3)
        e_type[i] = electronType()

        # Choose the z component of the momentum (step 4)
        if e_type[i] == 0: #'core':
            # assign energy based on gaussian distribution
            z_energy[i] = random.gauss(1, 4) # mu=1, sigma=4 eV

        elif e_type[i] == 1: #'valence':
            z_energy[i] = np.nan
            # assign energy based on fermi momentum distribution
            # TODO

        else:
            print('Something went very, very wrong')

        # Calculate the angle??? (step 7 in pseudo-pseudo-code in OneNote)
        # math in the OneNote, Tuesday Feb 8th notes
        theta[i] = m.atan(z_energy[i] / PHOTON_ENERGY ) # 511 keV, guassian in eV

        # find height that particle will hit at
        height[i] = m.tan(theta[i]) * DETECTOR_DIST

    ########################################################
    # Plot stuff
    ########################################################

    # plot the gaussian dist generated
    # plt.figure()
    # plt.hist(z_energy)
    # plt.show()

    # plot the fermi dist generated

    # histogram height on x axis
    plt.figure()
    plt.hist(height, bins=numbins_height)
    plt.show()
    
    # histogram for each electron type
    






if __name__=="__main__":
    main()