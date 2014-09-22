#!/usr/bin/python
#include [[functions_smsc.py]]
import functions_smsc as of 
import logging ,sys, getopt, OPA
Hartree2cm_1=219474.63 

def usage():
   print("usage: smallscipt.py[-o spectfile | -l logging | -t Temp | -e Energy]  arg1 arg2 ")
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
   """ script for the spectrum-calculation using two log-files from (see usage)
   g09 frequency calculations (input-argument 1 (initial state) and 2(final state))

   The program requires the non-standard libraries
   numpy
   matplotlib (this can be excluded easily, see README)
   """

   #===handling of input-arguments:===
   #Input(sys.argv[2]) #just ignore further argements!
   spectfile='/dev/null' # by default: discart broadened spectrum
   E0=0
   T=300
   if argv is None:
      argv = sys.argv[1:]
   try:
      opts, args=getopt.getopt(argv, 'h:o:l:t:e:', ["help", "logging=" ,"out=", "Temperature=", "Energy="])
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
      print "log-files",args," are miss-typed or don't exist.\n"
      print opts
      usage()
      sys.exit(2)
   for opt,s in opts:
      if opt in ['-o', '--out']:
	 spectfile=s
      elif opt in ['-l', '--logging']:
	 invokeLogging(s)
      elif opt in ['-t','--Temperature']:
	 T=float(s)
      elif opt in ['-e','--Energy']:
	 E0=float(s) #test, whether s is an int
   if ['-l', '--logging'] not in opts:
      invokeLogging()
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
   #Gauf=of.gaussianfreq(ContntInfo, dim) #read frequencies calculated by g09 from file
   #Gauf/=Hartree2cm_1  #convert to atomic units
   #F, P, CartCoord=of.TrafoCoord(F, P, CartCoord, dim)
   logging.info('Cartesion coordinates of final system:\n'+repr(CartCoord[0].T)+'\n'+ repr(CartCoord[1].T))
   logging.info('forces:\n'+repr(F[0])+'final state:\n'+repr(F[0]))
   
   #=== Calculate Frequencies and normal modes ===
   L, f, Lsorted=of.GetL(dim, mass,F, P)
   #L2=of.extractL(ContntInfo, dim)
   #of.replace(inputs[0], f[1], Lsorted[1])
   
   ##=== Spectrum calculation ===
   #J, K=of.Duschinsky(L2, mass, dim, CartCoord) #use gaussians normal modes
   J, K=of.Duschinsky(Lsorted, mass, dim, CartCoord) #use own quantities
   print 'Energies:', Energy[0], Energy[1]
   print '0-0-transition:', (Energy[0]-Energy[1])*Hartree2cm_1

   ##==calculate HR-Spectrum==
   HR, funi= of.HuangR(K, f)
   #linspect=of.calcspect(HR, funi, Energy[0]-Energy[1], E0, 5, 5, T, "TPA")
   linspect=of.calcspect(HR, funi, Energy[0]-Energy[1], E0, 5, 5, T)
   #of.outspect(61009, linspect, 3, spectfile)

   #==calculate Duschinsky-Rotated Spectrum taking OPA into account==
   #linspect=OPA.FCfOPA(J,K,f,Energy[0]-Energy[1],4, T, E0)
   #of.outspect(61009, linspect, 3, spectfile)
   
if __name__ == "__main__":
   main(sys.argv[1:])
   
