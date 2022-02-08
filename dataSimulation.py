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

def electronType() -> str:
    return random.choice([0,1]) # 0=core, 1=valence

def main():

    e_type = np.empty(N)        # electron type
    z_momentum = np.empty(N)    # momentum in z direction

    for i in range(0, N):
        # Get random choice of valence or core electron hit by positron (step 3)
        e_type[i] = electronType()

        # Choose the z component of the momentum (step 4)
        if e_type[i] == 0: #'core':
            # assign momentum based on gaussian distribution
            z_momentum[i] = random.gauss(1, 4) # mu=1, sigma=4 eV

        elif e_type[i] == 1: #'valence':
            z_momentum[i] = np.nan
            # assign momentum based on fermi momentum distribution
            # TODO

        else:
            print('Something went very, very wrong')

        # Calculate the angle??? (step 7 in pseudo-pseudo-code in OneNote)





    ########################################################
    # Plot stuff
    ########################################################

    # plot the gaussian dist generated
    plt.figure()
    plt.hist(z_momentum)
    plt.show()

    # plot the fermi dist generated
    
    






if __name__=="__main__":
    main()