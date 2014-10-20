#!/usr/bin/python
# filename: broadening.py
import numpy as np, re, logging
# Below are the conversion factors and fundamental constant

def handel_input(opt):
   gridfile=None
   gamma=1 #by default: only slight broadening
   gridpt=5000
   omega=None
   minfreq=0
   maxfreq=0
   shape='g'
   spectfile=None

   tmpgrid=re.findall(r"(?<=grid=)[ \=\s\w\.,]+", opt, re.M)
   if len(tmpgrid)==1: 
      # i.e. if grid is specified
      grid=re.findall(r"[\w\.]+", tmpgrid[0], re.M)
      if len(grid)==1:
	 try:
	    gridpt=float(grid[0])
	 except ValueError: # if grid is no a number
	    gridfile=grid[0]
      elif len(grid)==3:
	 gridpt=float(grid[0])
	 minfreq=float(grid[1])
	 maxfreq=float(grid[2])
      if gridfile!=None:
	 #read file in format of linspect
	 grid=[]
	 with open(gridfile) as f:
	    lis=[line.split() for line in f]  # create a list of lists
	    for i,x in enumerate(lis):        #print the list items 
	       grid.append(float(x[0]))
	 omega=np.zeros(len(grid))
	 for i in range(len(grid)):
	    omega[i]=grid[i]
   if (re.search(r"(?<=gamma=)[ \d\.,]+", opt, re.I) is not None)  is True:
      gamma=float(gamma[0])
   shape=re.findall(r"(?<=shape=)[ \w]+", opt, re.I)
   if len(shape)>0:
      if shape[0] in ["lorentzian", "Lorentzian", "L", "l"]:
	 shape="l"
      elif shape[0] in ["gaussian", "Gaussian", "G", "g"]:
	 shape="g"
   return omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape

def OPA2TPA(OPAfreq,freq00, OPAintens,intens00, mode):
   length=len(OPAfreq)
   TPAfreq=np.zeros((length+1)*(length+2)//2+length+1) #this is overestimation of size...
   TPAintens=np.zeros((length+1)*(length+2)//2+length+1)
   TPAintens[0]=intens00 #this is OPA-part
   TPAfreq[0]=freq00
   #print intens00, freq00 , 0
   ind=1
   for i in range(length):
      TPAintens[ind]=OPAintens[i] #this is OPA-part
      TPAfreq[ind]=OPAfreq[i]
      #print TPAintens[ind], TPAfreq[ind], ind
      ind+=1
      for j in range(i+1,length):
	 if mode[i]==mode[j]: #both have same mode...
	    continue
	 TPAintens[ind]=OPAintens[i]*OPAintens[j]/intens00
	 TPAfreq[ind]=OPAfreq[i]+OPAfreq[j]-freq00
	 #print TPAintens[ind], TPAfreq[ind], ind
	 ind+=1
   return TPAfreq, TPAintens

def OPA23PA(OPAfreq,freq00, OPAintens,intens00, mode):
   length=len(OPAfreq)
   TPAfreq=[]#np.zeros((((length+1)*(length+2)+3)*(length+3))//6+length+1) #this is overestimation of size...
   TPAintens=[]
   TPAintens.append(intens00) #this is OPA-part
   TPAfreq.append(freq00)
   #print intens00, freq00 , 2
   for i in range(length):
      TPAintens.append(OPAintens[i]) #this is OPA-part
      TPAfreq.append(OPAfreq[i])
      #print TPAintens[-1], TPAfreq[-1] ,2
      for j in range(i+1,length):
	 if mode[i]==mode[j]: #both have same mode...
	    continue
	 TPAintens.append(OPAintens[i]*OPAintens[j]/intens00)
	 TPAfreq.append(OPAfreq[i]+OPAfreq[j]-freq00)
	 #print TPAintens[-1], TPAfreq[-1], 2
	 for k in range(j+1,length):
	    if mode[k]==mode[j]: #both have same mode...
	       continue
	    if mode[k]==mode[i]:
	       continue
	    TPAintens.append(OPAintens[i]*OPAintens[j]*OPAintens[k]/(intens00*intens00))
	    TPAfreq.append(OPAfreq[i]+OPAfreq[k]+OPAfreq[j]-2*freq00)
	    #print TPAintens[-1], TPAfreq[-1], 2
   freq=np.zeros(len(TPAfreq))
   intens=np.zeros(len(TPAintens))
   for i in range(len(freq)): #this can not be done by np.matrix() due to dimensionality...
	freq[i]=TPAfreq[i]
	intens[i]=TPAintens[i]
   return freq, intens

def outspect(T, opt, linspect, E=0):
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
   omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape =handel_input(opt)
   #read file in format of linspect

   #sort spectrum with respect to size of elements
   index=np.argsort(linspect[1],kind='heapsort')
   linspect[1]=linspect[1][index]
   linspect[2]=linspect[2][index]
   linspect[0]=linspect[0][index]
   #find transition with minimum intensity to be respected
   minint=0
   for i in range(len(linspect[1])):
      if linspect[1][i]>=0.001*linspect[1][-1]:
	 minint=i
	 break
   logging.warning('neglect '+repr(minint)+' transitions, use only '+repr(len(linspect[1])-minint)+" instead.")
   logging.info('minimal and maximal intensities:\n'+repr(linspect[1][minint])+' '+repr(linspect[1][-1]))

   #make TPA from OPA:
   if (re.search(r"to ?PA", opt, re.I) is not None) is True:
      n=re.findall(r"(?<=(to )[\d](?=PA)", opt, re.I)
      if n[0]=='2':
	 TPAfreq, TPAintens=OPA2TPA(linspect[0][minint:],E , linspect[1][minint:], 10, linspect[2][minint:])
	 minint=0
	 for i in range(len(TPAintens)):
	    if TPAintens[i]>=0.0001*TPAintens[-1]:
	       minint=i
	       break
	 logging.warning('for TPA: again neglect '+repr(minint)+
		     ' transitions, use only '+repr(len(TPAintens)-minint)+" instead.")
	 index=np.argsort(TPAintens,kind='heapsort')
	 TPAintens=TPAintens[index] #resort by intensity
	 TPAfreq=TPAfreq[index]
      elif n[0]=='3':
	 TPAfreq, TPAintens=OPA23PA(linspect[0][minint:],E , linspect[1][minint:], 10, linspect[2][minint:])
	 minint=0
	 for i in range(len(TPAintens)):
	    if TPAintens[i]>=0.0001*TPAintens[-1]:
	       minint=i
	       break
	 logging.warning('for 3PA: again neglect '+repr(minint)+
		     ' transitions, use only '+repr(len(TPAintens)-minint)+" instead.")
	 index=np.argsort(TPAintens,kind='heapsort')
	 TPAintens=TPAintens[index] #resort by intensity
	 TPAfreq=TPAfreq[index]
      elif n[0]==1:
	 TPAfreq=linspect[0][minint:]
	 TPAintens=linspect[1][minint:]
      else:
	 logging.critical("to <n>PA was given but not recognised.")
   else:
      TPAfreq=linspect[0][minint:]
      TPAintens=linspect[1][minint:]

   #find transition with minimum intensity to be respected
   logging.info("intensity, frequency,   2")
   for i in range(len(TPAfreq)):
	logging.info(repr(TPAintens[i])+"  "+repr(TPAfreq[i])+"  "+repr(2))

   #the range of frequency ( should be greater than the transition-frequencies)
   if omega==None:
      if minfreq==0:
	 minfreq=np.min(TPAfreq)-20-gamma
      if maxfreq==0:
	 maxfreq=np.max(TPAfreq)+20+gamma
   else:
      minfreq=omega[0]
      maxfreq=omega[-1]
   logging.warning('maximal and minimal frequencies: '+repr(maxfreq)+"  "+repr(minfreq))
   #truncate arrays and sort by index for further efficient processes
   #if no other grid is defined: use linspace in range
   if omega==None:
      omega=np.linspace(minfreq,maxfreq,gridpt)
      logging.info("omega is equally spaced")
   spect=np.zeros(len(omega))
   sigma=gamma*2/2.355 #if gaussian used: same FWHM

   index=np.argsort(TPAfreq,kind='heapsort') #sort by freq
   freq=TPAfreq[index]
   intens=TPAintens[index]
   mini=0
   for i in range(1,len(freq)):
      if freq[i]>=5*gamma+freq[0]:
	 maxi=i
	 break
   if spectfile==None:
      spectfile="/dev/null" #discart spectrum
   out = open(spectfile, "w")
   log = open("calculation.log", "a")
   log.write("broadened spectrum:\n frequency      intensity\n")
   #logging.critical("broadened spectrum:\n frequency      intensity")
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
	 log.write(u" {0}  {1}\n".format(omega[i] ,spect[i]))
	 #logging.critical(u" {0}  {1}".format(omega[i] ,spect[i]))
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
	 log.write(u" {0}  {1}\n".format(omega[i] ,spect[i]))
	 #logging.critical(u" {0}  {1}".format(omega[i] ,spect[i]))
   log.close()
   out.close()

version=1.1
# End of broadening.py
