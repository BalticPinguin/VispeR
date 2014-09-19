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
   gridpt=61009
   omega=None
   minfreq=0
   maxfreq=0
   shape='g'
   E00=0

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
      elif opt in ['-E','--energy']:
	 E00=float(s)
      elif opt in ['-s','--shape']:
	 if s in ['gaussian', 'g' ,'gauss', 'Gauss', 'Gaussian']:
	    shape='g'
	 elif s in ['lorentzian', 'l' ,'lorentz', 'Lorentz', 'Lorentzian']:
	    shape='l'
	 else:
	    print 'shape unknown. Please gaussian or lorentzian instead!'
   #else: #if no input is defined: senseless
   #   usage()
   #   sys.exit(2)
   if linspectrum==None:
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
   return linspectrum, omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape, E00

def OPA2TPA(OPAfreq,freq00, OPAintens,intens00, mode):
   length=len(OPAfreq)
   TPAfreq=np.zeros((length+1)*(length+2)//2+length+1) #this is overestimation of size...
   TPAintens=np.zeros((length+1)*(length+2)//2+length+1)
   TPAintens[0]=intens00 #this is OPA-part
   TPAfreq[0]=freq00
   print intens00, freq00 
   ind=1
   for i in range(length):
      TPAintens[ind]=OPAintens[i] #this is OPA-part
      TPAfreq[ind]=OPAfreq[i]
      print TPAintens[ind], TPAfreq[ind]
      ind+=1
      for j in range(i+1,length):
	 if mode[i]==mode[j]: #both have same mode...
	    continue
	 TPAintens[ind]=OPAintens[i]*OPAintens[j]/intens00
	 TPAfreq[ind]=OPAfreq[i]+OPAfreq[j]-freq00
	 print TPAintens[ind], TPAfreq[ind]
	 ind+=1
   return TPAfreq, TPAintens

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
      opts, args=getopt.getopt(argv, 'h:o:g:i:p:m:M:f:s:E:', 
           ["help","out=", "gamma=", "input=", "gridpt=", "minfreq=", "maxfreq=", "gridfile=","shape=","energy="])
   except getopt.GetoptError as err:
      print(err)
      usage()
      sys.exit(2)
   linspectrum, omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape, E00=handel_input(opts, args)
      
   #read file in format of linspect
   freq=[]
   intens=[]
   mode=[]
   
   with open(linspectrum) as f:
      lis=[line.split() for line in f]  # create a list of lists
      for i,x in enumerate(lis):        #print the list items 
         freq.append(float(x[0]))
         intens.append(float(x[1]))
	 mode.append(float(x[2]))
   linspect=np.zeros((3,len(freq)))
   linspect[0]=np.matrix(freq)
   linspect[1]=np.matrix(intens)
   linspect[2]=np.matrix(mode)

   #sort spectrum with respect to size of elements
   index=sort(linspect[1])
   linspect[1]=linspect[1][index]
   linspect[2]=linspect[2][index]
   linspect[0]=linspect[0][index]
   #find transition with minimum intensity to be respected
   minint=0
   for i in range(len(linspect[1])):
      if linspect[1][i]>=0.0001*linspect[1][-1]:
	 minint=i
	 break
   print 'neglect',minint,'transitions, use only', len(linspect[1])-minint, "instead."
   print('minimal and maximal intensities:\n', linspect[1][minint], linspect[1][-1])

   TPAfreq=linspect[0][minint:]
   TPAintens=linspect[1][minint:]
   ##make TPA from OPA:
   #TPAfreq, TPAintens=OPA2TPA(linspect[0][minint:], 18311.8877, linspect[1][minint:], 100, linspect[2][minint:])
   #index=sort(TPAintens)
   #TPAintens=TPAintens[index] #resort by intensity
   #TPAfreq=TPAfreq[index]

   #find transition with minimum intensity to be respected
   minint=0
   for i in range(len(TPAintens)):
      if TPAintens[i]>=0.0001*TPAintens[-1]:
	 minint=i
	 break
   print 'for TPA: again neglect',minint,'transitions, use only', len(TPAintens)-minint, "instead."

   #the range of frequency ( should be greater than the transition-frequencies)
   if omega==None:
      if minfreq==0:
	 minfreq=np.min(TPAfreq[minint:])-20-gamma
      if maxfreq==0:
	 maxfreq=np.max(TPAfreq[minint:]) +20+gamma
   else:
      minfreq=omega[0]
      maxfreq=omega[-1]
   print('maximal and minimal frequencies:\n', maxfreq, minfreq)
   #truncate arrays and sort by index for further efficient processes
   #if no other grid is defined: use linspace in range
   if omega==None: 
      omega=np.linspace(minfreq,maxfreq,gridpt)
      print "omega is linspaced"
   spect=np.zeros(len(omega))
   sigma=gamma*2/2.355 #if gaussian used: same FWHM

   index=sort(TPAfreq) #sort by freq
   freq=TPAfreq[index]
   intens=TPAintens[index]
   mini=0
   for i in range(1,len(freq)):
      if freq[i]>=5*gamma+freq[0]:
	 maxi=i
	 break
   out = open(spectfile, "w")
   if shape=='g':
      for i in range(len(omega)): 
     	 for j in range(maxi,len(freq)):
   	    if freq[j]>=5*gamma+omega[i]:
   	       maxi=j
   	       break
	 for j in range(max(0,mini),maxi):
   	    if freq[j]>=omega[i]-5*gamma:
   	       mini=j-1
   	       break
	 spect[i]=sum(intens[j]/(np.sqrt(2*np.pi)*sigma)*\
		  np.exp(-(omega[i]-freq[j])*(omega[i]-freq[j])/(2*sigma*sigma))
		  for j in range(mini, maxi))
	 out.write(u" {0}  {1}\n".format(omega[i] ,spect[i]))
   else:  #shape=='l':
      for i in range(len(omega)): 
	 for j in range(maxi,len(freq)):
   	    if freq[j]>=5*gamma+omega[i]:
   	       maxi=j
   	       break
     	 for j in range(max(0,mini),maxi):
   	    if freq[j]>=omega[i]-5*gamma:
   	       mini=j-1
   	       break
	 spect[i]=sum(intens[j]/np.pi*gamma/((omega[i]-freq[j])*(omega[i]-freq[j])+gamma*gamma)
		  for j in range(mini, maxi))
	 out.write(u" {0}  {1}\n".format(omega[i] ,spect[i]))

   out.close()

if __name__ == "__main__":
   outspect(sys.argv[1:])

version=0.0
# End of broadening.py
