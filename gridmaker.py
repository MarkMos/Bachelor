## This script creates 3-dimensional grids of k-vectors, delta1 and theta1 from
## k=-20 to 20 Mpc^-1 in all three dimensions and delta2 and theta2 from k=-20
## to 20 Mpc^-1 in all three dimensions. The spacing between grid points is 0.5 Mpc^-1
## in all three dimensions. It also produces a sum of deltas 1 and 2 and thetas
## 1 and 2 in this grid.
##
## Other than standard Python packages, this script uses two small modules I have
## written myself, FF2 and myfuncs as well as the classy wrapper for CLASS. CLASS
## is available at class-code.net

#Import necessary packages
import numpy as np
from classy import Class
from scipy.interpolate import interp1d
import FF2
import myfuncs as my
from numpy import linalg as LA
import pickle

#Run ClASS-code and get transfer function for delta_cdm
cosmo = Class()
cosmo.set({'output':'tCl mPk dTk vTk','P_k_max_1/Mpc':35, 'gauge':'Synchronous','z_pk':'100., 0.'})
cosmo.compute()
tr = cosmo.get_transfer(49)
delta_cdm = tr['d_cdm']
k= tr['k (h/Mpc)']*cosmo.h()

# Create interpolating function for transfer function
dfunc = interp1d(k, delta_cdm, kind='cubic', bounds_error=False, fill_value=0)

# define class for saving data
class data:

    def __init__(self, Dat):
        self.Dat = Dat
        self.K = []
        self.Kused = []
        self.Kl = []
        self.Klu = []
        self.d1 = []
        self.d2 = []
        self.t1 = []
        self.t2 = []
        self.dsum = []
        self.tsum = []

#define dataset
dataset = data('DatSetName')

#define array of k-vectors
KARRAY = my.grid3D(-20,20,81)
Kused = KARRAY[20:61,20:61,20:61,:]

dataset.K = KARRAY
dataset.Kused = Kused

#create array of random numbers fulfilling symmetry requirements
R = np.random.normal(0,1,(81,81,81))+1j*np.random.normal(0,1,(81,81,81))

for x in range(0,81):
    for y in range(0,81):
        for z in range(0,81):
            R[x,y,z] = np.conj(R[80-x,80-y,80-z])

Rused = R[20:61,20:61,20:61]

#creating array with lengths of k-vectors, setting length-0 to 1e-9 to prevent division with zero
Klen = LA.norm(KARRAY, axis=3)
Klen[40,40,40] = 1e-9
Klu = Klen[20:61,20:61,20:61]

#calculate physical delta1 and theta1 from CLASS transfer function
As=np.exp(3.089)*1e-10
ns=0.9655
d1 = dfunc(Klen)*R*np.sqrt(As/(2*np.pi**2))*Klen**(-3./2.) * (Klen/0.05)**((ns-1)/2)
d1used = d1[20:61,20:61,20:61]
theta1 = -d1
theta1used = -d1used

dataset.d1 = d1
dataset.t1 = theta1

#creating arrays for delta2 and theta2
d2 = np.empty((41,41,41), dtype=complex)
theta2 = np.empty((41,41,41), dtype=complex)

#calculating delta and theta2
d1mK = np.flip(np.flip(np.flip(d1,2) ,1) ,0)
for x in range(0,41):
    for y in range(0,41):
        for z in range(0,41):
            K = Kused[x,y,z,:]
            KminQ = K - Kused
            KmQl = LA.norm(KminQ, axis=3)
            KmQl[KmQl==0]=1e-9
            cPhi = np.sum(KminQ*Kused,3)/(Klu*KmQl)
            F2 = FF2.F2s(Klu,KmQl,cPhi)
            G2 = FF2.G2s(Klu,KmQl,cPhi)
            intfuncd = F2*d1used*d1mK[40-x:81-x,40-y:81-y,40-z:81-z]
            intfunct = G2*d1used*d1mK[40-x:81-x,40-y:81-y,40-z:81-z]
            d2[x,y,z] = np.sum(intfuncd)*0.5**3
            theta2[x,y,z] = -np.sum(intfunct)*0.5**3

dataset.d2 = d2
dataset.t2 = theta2
#creating full deltas and thetas
d1pd2 = d1[20:61,20:61,20:61]+d2
t1pt2 = theta1[20:61,20:61,20:61]+theta2

dataset.dsum = d1pd2
dataset.tsum = t1pt2

#saving data with pickle
pickle.dump(dataset,open("data.p","wb"))
