#!/usr/bin/python
# filename: broadening.py
import numpy as np, getopt, sys
from copy import deepcopy # for function sort(). Probably find a better function!!
# Below are the conversion factors and fundamental constant

def usage():
   print "usage: broadening.py [-o file to write in | -g broadening-parameter | -i input-file | -p number of gridpoints",\
	    "| -m minimal frequency | -M maximal frequency]"
   print "alternative usage: broadening.py [-o file to write in | -g broadening-parameter | -i input-file | -gf gridfile]"
   print "if no out- or input-file is specified, the program exits."

def sort(f):
   """
   This function sorts float-arrays by absolute value (lowest argument first arguments)

   **PARAMETERS:**
   f:  array to be sorted

   **RETURNS:**
   order of indices of f

   **NOTE:**
   the elements of f should not exceed 3e300 (otherwise the sorting will fail) 
   of indices sorted by the size of respective elements.
   The sorting by the size of elemens in A (largest first) can be done by 

   index=sort(A) 
   B=A[index]

   where B will be the sorted array.
   """
   index=np.zeros(len(f), dtype=int)
   tmp=deepcopy(f) #for avoiding side effects
   for i in range(len(f)):
      index[i]=np.argmin(np.abs(tmp)) # there can be frequencies < 0 as well...
      tmp[index[i]]=5e+300 # this can be considered as smaller than all elements...
   return index

def handel_input(opts, args):
   spectfile=None
   gridfile=None
   linspectrum=None
   gamma=1 #by default: only slight broadening
   gridpt=6000
   omega=None
   minfreq=0
   maxfreq=0

   if opts in ['-h', '--help']:
      usage()
      sys.exit()
   for opt,s in opts:
      print opt,s
      if opt in ['-o', '--out']:
	 spectfile=s
      elif opt in ['-g','--gamma']:
	 gamma=float(s)
      elif opt in ['-i','--input']:
	 linspectrum=s 
      elif opt in ['-f','--gridfile']:
	 gridfile=s 
      elif opt in ['-p','--gridpt']:
	 gridpt=float(s)
      elif opt in ['-m','--minfreq']:
	 minfreq=float(s)
      elif opt in ['-M','--maxfreq']:
	 maxfreq=float(s)
   if spectfile!=None:
      out = open(spectfile, "w")
   else: #if no input is defined: senseless
      usage()
      sys.exit(2)
   if linspectrum!=None:
      #read file in format of linspect
      freq=[]
      intens=[]
      with open(linspectrum) as f:
	 lis=[line.split() for line in f]  # create a list of lists
	 for i,x in enumerate(lis):        #print the list items 
	    freq.append(float(x[0]))
	    intens.append(float(x[1]))
      linspect=np.zeros((2,len(freq)))
      linspect[0]=np.matrix(freq)
      linspect[1]=np.matrix(intens)
   else: #if no input is defined: senseless
      usage()
      sys.exit(2)
   if gridfile!=None:
      #read file in format of linspect
      grid=[]
      with open(gridfile) as f:
	 lis=[line.split() for line in f]  # create a list of lists
	 for i,x in enumerate(lis):        #print the list items 
	    grid.append(float(x[0]))
      omega=np.zeros(len(grid))
      print np.shape(omega)
      for i in range(len(grid)):
	 omega[i]=grid[i]
      print np.shape(omega)
   return linspect, omega, out, gamma, gridpt, minfreq, maxfreq

def outspect(argv):
   """This function calculates the broadened spectrum given the line spectrum, 
   frequency-rage and output-file whose name is first argument. 
   As basis-function a Lorentzian is assumed with a common width.
   
   **PARAMETERS:**
   spectfile: file, the result is written in (ascii-table). 
   In addition a graph is created and shown on the fly. This graph is not saved.
   gridpt:    number of grid-points to be used for the calculation
   linspect:  line-spectrum list (frequency, intensity) 
   gamma:     broadening constant for the Lorentzians. It is the same for all peaks
   
   """
   if argv is None:
      argv = sys.argv[1:]
   try:
      opts, args=getopt.getopt(argv, 'h:o:g:i:p:m:M:f:', ["help","out=", "gamma=", "input=", "gridpt=", "minfreq", "maxfreq", "gridfile"])
   except getopt.GetoptError as err:
      print(err)
      usage()
      sys.exit(2)
   linspect, omega, out, gamma,gridpt,minfreq,maxfreq=handel_input(opts, args)
   #sort spectrum with respect to size of elements
   index=sort(linspect[1])
   linspect[1]=linspect[1][index]
   linspect[0]=linspect[0][index]
   #find transition with minimum intensity to be respected
   minint=0
   for i in range(len(linspect[1])):
      if linspect[1][i]>=0.001*linspect[1][-1]:
	 minint=i
	 break
   print 'neglect',minint,'transitions, use only', len(linspect[1])-minint, "instead."
   print('minimal and maximal intensities:\n', linspect[1][minint], linspect[1][-1])
   #the range should be greater than the transition-frequencies
   if omega==None:
      if minfreq==0:
	 minfreq=np.min(linspect[0][minint:])-20-gamma
	 minfreq-=20+gamma
      if maxfreq==0:
	 maxfreq=np.max(linspect[0][minint:]) +20+gamma
   else:
      minfreq=omega[0]
      maxfreq=omega[-1]
   print('maximal and minimal frequencies:\n', maxfreq, minfreq)
   #if no other grid is defined: use linspace in range
   if omega==None: 
      omega=np.linspace(minfreq,maxfreq,gridpt)
      print "omega is linspaced"
   spect=np.zeros(len(omega))
   #only those with high-enough intensities are respected
   for i in range(len(omega)): 
      spect[i]=sum(linspect[1][j]/np.pi*gamma/((omega[i]-linspect[0][j])*(omega[i]-linspect[0][j])+gamma*gamma)
	    for j in range(minint,len(linspect[1])))
      print i, omega[i],spect[i]
      out.write(u" {0}  {1}\n".format(omega[i] ,spect[i]))

if __name__ == "__main__":
   outspect(sys.argv[1:])

version=0.0
# End of broadeninge.py
