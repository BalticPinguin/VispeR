#!/usr/bin/python
#include [[functions_smsc.py]]
import functions_smsc as of 
import logging ,sys

""" script for the spectrum-calculation using two log-files from
g09 frequency calculations (input-argument 1 (initial state) and 2(final state))
"""
logging.basicConfig(filename='calculation.log',level=logging.DEBUG)
logging.info('Initializing the log-file')

#===handling of input-arguments:===
inputs =[]
inputs.append(sys.argv[1])
inputs.append(sys.argv[2])
logging.debug('indut-data: (number of files, names) '+repr(len(inputs))+' '+repr(inputs))
spectfile=sys.argv[3] #this should be a non-abligatory argument
logging.info('START of calculations. initial-state file: '+\
	 repr(inputs[0])+', final-state file: '+repr(inputs[1])+\
	 '. The spectrum will be printed into '+repr(spectfile)+'.')

# ==look for the investigated molecule and where opt/freq was researched ==
ContntInfo=of.Contentcheck(inputs) # function: makes tests referring integrity, gathers general information
if len(ContntInfo)!=2:
   logging.error('one of the files has invalid data.')
assert len(ContntInfo)==2, 'one of the files has invalid data.'
logging.debug("Content info:\n"+ repr(ContntInfo))

for i in range(2): #read coordinates, force constant, binding energies from log-files and calculate needed quantities
   dim, Coord, mass, B, A, E=of.ReadLog(ContntInfo[i][0])
   if i is 0:# do only in first run
      F, CartCoord, X, P, Energy=of.quantity(dim)
   logging.debug("Dimensions: "+ repr(dim)+ '\n Masses: '+ repr(mass**2))
   X[i]=B
   F[i]=A
   Energy[i]=E
   CartCoord[i]=Coord
   P[i]=of.GetProjector(X[i], dim, mass, CartCoord[i])
   logging.debug('Projector onto internal coordinate subspace\n'+ repr(P[i]))
#print 'Force-difference\n',F[0]-F[1]
logging.info('difference of minimum energy between states:  '+ repr(Energy[1]-Energy[0]))
Gauf=of.gaussianfreq(ContntInfo, dim) #read frequencies calculated by g09 from file
Gauf/=219474.63  #convert to atomic units

#F, P, CartCoord=of.TrafoCoord(F, P, CartCoord, dim)
logging.info('Cartesion coordinates of final system:\n'+repr(CartCoord[0].T)+'\n'+ repr(CartCoord[1].T))
logging.info('forces:\n'+repr(F[0])+'final state:\n'+repr(F[0]))

N, L, f, Lsorted=of.GetL(dim, mass,F, P)
#L1, f1, Lsorted1=of.GetL1(dim, mass, F, Gauf, P)
L2=of.extractL(ContntInfo, dim)

#=== Spectrum calculation ===
J, K=of.Duschinsky(N, L, mass, dim, CartCoord)
#==calculate HR-Spectrum==
HR_unif, funi, HR_multif, sortfG, sortfE= of.HuangR(K, f)
linspect=of.calcspect(HR_unif[-5:], funi[-5:], Energy[1]-Energy[0], 5, 5)
#of.outspect(spectfile, 3000, linspect, 4002)

#==calculate Duschinsky-Rotated Spectrum==
linespect=of.FCf(J,K,f,Energy[1]-Energy[0],5)
#of.outspect(spectfile, 3000, linspect, 80)

logging.info('END of calculations')

