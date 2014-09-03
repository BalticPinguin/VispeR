#!/usr/bin/python
#include [[functions_smsc.py]]
import functions_smsc as of 
import logging ,sys, getopt
Hartree2cm_1=219474.63 

def usage():
   print("usage: smallscipt.py[-o spectfile | -l logging | -t Temp]  arg1 arg2 ")
   print("spectfile: (not nonobligatory) filename for spectrum-information")
   print("logging:   (default: error) debug-mode")
   print("Temp:      (default: 300K) at the moment not used")
   print("arg1:      initial-state file")
   print("arg2:      final-state file\n")
   print("or just call with '-h' to see this help")

def invokeLogging(mode="ERROR"):
   if mode in ['debug','DEBUG','Debug']:
      logging.basicConfig(filename='calculation.log',level=logging.DEBUG)
   elif mode in ['info','INFO','Info']:
      logging.basicConfig(filename='calculation.log',level=logging.INFO)
   elif mode in ['error', 'ERROR','Error']:
      logging.basicConfig(filename='calculation.log',level=logging.ERROR)
   else:
      logging.basicConfig(filename='calculation.log',level=logging.ERROR)
      logging.error("logging-mode not recognized. Using 'error' instead")
   logging.info('Initializing the log-file')

def main(argv=None):
   """ script for the spectrum-calculation using two log-files from
   g09 frequency calculations (input-argument 1 (initial state) and 2(final state))

   The program requires the non-standard libraries
   numpy
   matplotlib (this can be excluded easily)
   """

   #===handling of input-arguments:===
   if argv is None:
      argv = sys.argv[1:]
   try:
      opts, args=getopt.getopt(argv, 'h:o:l:t:', ["help", "logging=" ,"out=", "Temperature="])
   except getopt.GetoptError as err:
      print(err)
      usage()
      sys.exit(2)
   if opts in ['-h', '--help']:
      usage()
      sys.exit()
   if len(args)==2:
      inputs=args
   else:
      print ("log-files",args," are miss-typed or don't exist.\n")
      usage()
      sys.exit(2)
   for opt,s in opts:
      if opt in ['-o', '--out']:
	 spectfile=s
      elif opt in ['-l', '--logging']:
	 invokeLogging(s)
      elif opt in ['-t','--Temperature']:
	 T=s #test, whether s is an int
   if ['-l', '--logging'] not in opts:
      invokeLogging()
   if ['-t', '--Temperature'] not in opts:
      T=300
   T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K

   logging.debug('indut-data: (number of files, names) '+repr(len(inputs))+' '+repr(inputs))
   logging.info('START of calculations. initial-state file: '+\
	 repr(inputs[0])+', final-state file: '+repr(inputs[1])+\
	 '. The spectrum will be printed into '+repr(spectfile)+'.')
   
   # ==look for the investigated molecule and where opt/freq was researched ==
   ContntInfo=of.Contentcheck(inputs) # function: makes tests referring integrity, gathers general information
   if len(ContntInfo)!=2:
      logging.error('one of the files has invalid data.')
   assert len(ContntInfo)==2, 'one of the files has invalid data.'
   logging.debug("Content info:\n"+ repr(ContntInfo))
   for i in range(2): 
      #read coordinates, force constant, binding energies from log-files and calculate needed quantities
      dim, Coord, mass, B, A, E=of.ReadLog(ContntInfo[i][0])
      if i is 0:# do only in first run
	 F, CartCoord, X, P, Energy=of.quantity(dim) #create respective quantities (empty)
      logging.debug("Dimensions: "+ repr(dim)+ '\n Masses: '+ repr(mass**2))
      X[i],F[i],Energy[i]=B, A, E
      CartCoord[i]=Coord
      P[i]=of.GetProjector(X[i], dim, mass, CartCoord[i])
      logging.debug('Projector onto internal coordinate subspace\n'+ repr(P[i]))
   logging.info('difference of minimum energy between states:  '+ repr(Energy[1]-Energy[0]))
   Gauf=of.gaussianfreq(ContntInfo, dim) #read frequencies calculated by g09 from file
   Gauf/=Hartree2cm_1  #convert to atomic units
   print('energy-diff: ', Energy[0]-Energy[1])
   #F, P, CartCoord=of.TrafoCoord(F, P, CartCoord, dim)
   logging.info('Cartesion coordinates of final system:\n'+repr(CartCoord[0].T)+'\n'+ repr(CartCoord[1].T))
   logging.info('forces:\n'+repr(F[0])+'final state:\n'+repr(F[0]))
   
   #=== Calculate Frequencies and normal modes ===
   N, L, f, Lsorted=of.GetL(dim, mass,F, P)
   #L1, f1, Lsorted1=of.GetLstab(dim, mass,F, Gauf, P) ##at the moment not working/no sensible functionality
   L2=of.extractL(ContntInfo, dim)
   
   ##=== Spectrum calculation ===
   #J, K=of.Duschinsky(L2, mass, dim, CartCoord) #use gaussians normal modes
   J, K=of.Duschinsky(Lsorted, mass, dim, CartCoord) #use own quantities
   
   ##==calculate HR-Spectrum==
   #HR_unif, funi, HR_multif, sortfG, sortfE= of.HuangR(K, f)
   HR_unif, funi= of.HuangR(K, f)
   linspect=of.calcspect(HR_unif, funi, Energy[0]-Energy[1], 5, 5, T)
   of.outspect(3000, linspect, 20, spectfile)
   
   #==calculate Duschinsky-Rotated Spectrum==
   #of.FCf(J,K,f,Energy[0]-Energy[1],5)
   #of.outspect(3000, linspect, 20,spectfile)
   
   #==calculate Duschinsky-Rotated Spectrum using files==
   #of.fileFCf(J,K,f,Energy[0]-Energy[1],5)
   #of.sortfile()
   #of.fileoutspect(3000, linspect, 20,spectfile)
   
   logging.info('END of calculations')

if __name__ == "__main__":
   main(sys.argv[1:])
   
