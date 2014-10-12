#!/usr/bin/python
#include [[functions_smsc.py]]
import functions_smsc as of 
import OPA 
import broadening as br
import sys, os, logging, re, mmap, numpy as np

def usage():
   print "usage: smallscript <input-file>"

def invokeLogging(mode="important"):
   logging.basicConfig(format='%(message)s')
   if mode in ['all', 'ALL', 'All']:
      logging.basicConfig(filename='calculation.log',level=logging.DEBUG)
   elif mode in ['detailed', 'DETAILED', 'Detailed']:
      logging.basicConfig(filename='calculation.log',level=logging.INFO)
   elif mode in ['medium', 'MEDIUM','Medium']:
      logging.basicConfig(filename='calculation.log',level=logging.WARNING)
   elif mode in ['important', 'IMPORTANT','Important']:
      logging.basicConfig(filename='calculation.log',level=logging.ERROR)
   elif mode in ['short', 'SHORT','Short']:
      logging.basicConfig(filename='calculation.log',level=logging.CRITICAL)
   else:
      logging.basicConfig(filename='calculation.log',level=logging.ERROR)
      logging.error("logging-mode not recognized. Using 'error' instead")
   logging.info('Initializing the log-file')

def main(argv=None):
   assert len(argv)==1, 'exactly one argument required.'
   try:
      infile=open(argv[0], "r")
      f=mmap.mmap(infile.fileno(), 0, prot=mmap.PROT_READ)
      infile.close()
   except IOError:
      print "file", inputf, "not found."
      usage()
      return 2
   opts=[]
   todo=0
   if (re.search(r"HR-fact",f, re.I) is not None) is True:
      todo+=1
   opts.append(re.findall(r"(?<=HR-fact)[\w\d.,=() -]",f,re.I))
   if (re.search(r"FC-spect",f, re.I) is not None) is True:
      if (re.search(r"HR-file: ",f, re.I) is not None) is True:
	 #calculation of HR-facts not neccecary
	 todo=2
      else: #if 
	 todo=3
   opts.append(re.findall(r"(?<=FC-spect)[\w\d.=,() -]+",f,re.I))
   if (re.search(r"Duschinsky-spect",f, re.I) is not None) is True:
      if todo==0:
	 todo=5
      else:
	 todo+=4
   opts.append(re.findall(r"(?<=Duschinsky-spect)[\w\d=.(), -]",f,re.I))
   if (re.search(r"Broadening",f, re.I) is not None) is True:
      todo+=8
   opts.append(re.findall(r"(?<=Broadening\()[\w\d.,()= -](?==\))",f,re.I))
   if todo>=16 or todo==0 or todo in [4,6,9,13]:
      print "options for calculation don't make sense. Please check the input-file!"
      return 0
   print opts

   ############delete () from opts
   if np.mod(todo,2)==1: 
      #calculation up to HR-facts needed (FC- or Duschinsky spect)
      if opts[0]!=[]:
	 opt=opts[0][0]
      elif opts[1]!=[]:
	 opt=opts[1][0]
      elif opts[2]!=[]:
	 opt=opts[2][0]
      else:
	  print 'why am I here?'
	  return 2
      #invoke logging: take options of HR-fact, if exists, else: from FC-spect, else: from Duschinsky-spect
      loglevel=re.findall(r"(?<=level=)[\w]",opt, re.I)
      if loglevel==[]:
	 invokeLogging()
      else:
	 invokeLogging(loglevel)
      initial=re.findall(r"(?<=initial: )[\w.]+",f, re.I)
      final=re.findall(r"(?<=final: )[\w.]+",f, re.I)
      assert len(initial)>0, 'no initial state found!'
      assert len(final)>0, 'no final state found!'
      for i in range(len(initial)):
	 assert os.path.isfile(initial[i]) and os.access(initial[i], os.R_OK),\
	       initial[i]+' is not a valid file name or not readable.'
      for i in range(len(final)):
	 assert os.path.isfile(final[i]) and os.access(final[i], os.R_OK),\
	       final[i]+' is not a valid file name or not readable.'
      for i in range(len(initial)):
	 #read coordinates, force constant, binding energies from log-files and calculate needed quantities
	 dim, Coord, mass, B, A, E=of.ReadLog(initial[i])
	 if i is 0:# do only in first run
	    F, CartCoord, X, P, Energy=of.quantity(dim, len(initial)+len(final)) #creates respective quantities (empty)
	    logging.debug("Dimensions: "+ repr(dim)+ '\n Masses: '+ repr(mass**2))
	 X[i],F[i],Energy[i]=B, A, E
	 CartCoord[i]=Coord
	 P[i]=of.GetProjector(X[i], dim, mass, CartCoord[i])
	 logging.info('Projector onto internal coordinate subspace\n'+ repr(P[i]))
      for i in range(len(final)):
	 dim, Coord, mass, B, A, E=of.ReadLog(final[i]) 
	 X[i],F[i],Energy[i]=B, A, E
	 CartCoord[i]=Coord
	 P[i]=of.GetProjector(X[i], dim, mass, CartCoord[i])
	 logging.info('Projector onto internal coordinate subspace\n'+ repr(P[i]))
      for i in range(len(initial)):
	 for j in range(len(final)):
	    logging.warning('difference of minimum energy between states: '+ repr(Energy[j+len(initial)]-Energy[i]))
      for i in range(len(initial)):
	 logging.debug('Cartesion coordinates of initial state '+repr(i)+':\n'+repr(CartCoord[i].T)+'\n')
      for j in range(len(initial),len(final)+len(initial)):
	 logging.debug('Cartesion coordinates of final state '+repr(j)+':\n'+repr(CartCoord[j].T)+'\n')
      logging.info('forces:')
      for i in range(len(initial)):
	 logging.info(repr(i)+'. initial state: \n'+repr(F[i]))
      for j in range(len(initial), len(final)+len(initial) ):
	 logging.info(repr(j-len(initial))+'. final state: \n'+repr(F[j]))
   
      #Calculate Frequencies and normal modes
      L, f, Lsorted=of.GetL(dim, mass,F, P)
      J, K=of.Duschinsky(Lsorted, mass, dim, CartCoord)
      #calculate HR-spect
      HR, funi= of.HuangR(K, f)
      if (re.search(r"makeLog", opt, re.I) is not None) is True:  
	 for i in range(len(initial)): #### this needs to be enhanced
	    of.replace(initial[i], f[i], Lsorted[i])

   if np.mod(todo,4)>=2:
      ###calculate FC-spect
      opt=opts[1]
      if (re.search(r"broaden",opt, re.I) is not None) is True and todo<8:
	 todo+=8

   if np.mod(todo,8)>=4:
      ###calculate Duschinsky-sppect
      opt=opts[2]
      if (re.search(r"broaden",opt, re.I) is not None) is True and todo<8: 
	 todo+=8
   if np.mod(todo,16)>=8:
      ###calculate Broadening
      opt=opts[3]

   T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K
   
if __name__ == "__main__":
   main(sys.argv[1:])
