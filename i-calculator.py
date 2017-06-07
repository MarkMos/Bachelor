# Import necessary packages
import numpy as np
import pickle
from numpy import linalg as LA

# Declare data-class as used in my data
class data:

    def __init__(self, Dat):
        self.Dat = Dat
        self.K = []
        self.Kl = []
        self.Kused = []
        self.Klu = []
        self.d1 = []
        self.d2 = []
        self.t1 = []
        self.t2 = []
        self.dsum = []
        self.tsum = []

#Load data
dataset = pickle.load(open("data.p", "rb"))

#Define used data

Kdm = dataset.Kused
delta = dataset.dsum
theta = dataset.tsum
# If one does not wish to compare the full delta and theta, but rather the individual components, then use either of
#delta = dataset.d1[20:61,20:61,20:61]
#delta=dataset.d2
#theta = dataset.t1[20:61,20:61,20:61]
#theta = dataset.t2

#thetas is alwasy the full theta, as this is the theta with which the sign is correlated. The nonlinear mode-coupling is only seen when using the sign of the full theta
thetas = dataset.tsum

# Arrays for the binning of i^theta_theta,delta (abbreviated ittd) and the number of datapoints in each bin are correlated
ittd = np.zeros((18,18))
nums = np.zeros((18,18))

# In order to keep track of how far along the calculations has proceeded (the calculation takes several hours, I used 17-30 hours depending on pc), I have printed progress regularly. Here I set the number of ittds calculated so far to 0
numdone = 0

# Here I calculate my ittds
for xdm in range(0,41):
    #Here I print how far along the calculation is
    print numdone/41.0**6
    for ydm in range(0,41):
        for zdm in range(0,41):
            #Here I calculate the length of my k''-vector, so that it is only done once for each k''
            KdmL = LA.norm(Kdm[xdm,ydm,zdm,:])

            #Here I set the amount of calculations done to include all the calculations that will be done for this k''
            numdone = numdone + 41**3
            for xm in range(0,41):
                for ym in range(0,41):
                    for zm in range(0,41):
                        # Here I define my k as k'+k''
                        K = Kdm[xdm,ydm,zdm,:] + Kdm[xm,ym,zm,:]

                        # I determine if my k-vector falls outside the grid where my delta and theta are defined, if yes, I skip this calculation, as I have no theta(K) to use
                        if any(K>10) or any(K<-10):
                            continue
                        #Here I calculate the length of my k-vecor
                        KL = np.sqrt(K[0]*K[0]+K[1]*K[1]+K[2]*K[2])

                        #and here the length of my k'-vector squared
                        KmLsq = Kdm[xm,ym,zm,0]*Kdm[xm,ym,zm,0]+Kdm[xm,ym,zm,1]*Kdm[xm,ym,zm,1]+Kdm[xm,ym,zm,2]*Kdm[xm,ym,zm,2]

                        #I check if the length of k' is zero, if it is, I set it to a very small number to avoid division by zero
                        if KmLsq == 0:
                            KmLsq = 1e-9

                        #Here I calculate i_theta,delta
                        idt = (K[0]*Kdm[xm,ym,zm,0]+K[1]*Kdm[xm,ym,zm,1]+K[2]*Kdm[xm,ym,zm,2])/KmLsq * theta[xm,ym,zm]* delta[xdm,ydm,zdm]


                        # Here I calculate the appropriate index from my k-vector in order to find theta(k)
                        Ki = (K + 10)*2
                        thetaK = thetas[int(Ki[0]),int(Ki[1]),int(Ki[2])]

                        #Here I create integer indexes for binning by length of k and k''
                        iKdm = int(KdmL//1)
                        iK = int(KL//1)

                        #Here I finally calculate ittd and add it to the array, as well as adding 1 to the datapoint couter for this bin
                        ittd[iKdm,iK] = ittd[iKdm,iK] + np.real(idt) * np.sign(np.real(thetaK)) + np.imag(idt) * np.sign(np.imag(thetaK))
                        nums[iKdm,iK] = nums[iKdm,iK] +1

# The ittd and datapoint-counter arrays are saved
np.save('ittd',ittd)
np.save('nums',nums)
