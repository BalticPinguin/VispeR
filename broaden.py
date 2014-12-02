#!/usr/bin/python
# filename: broadening.py
import numpy as np, re
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

def OPA2nPA(OPAfreq,freq00, OPAintens, intens00, mode, n):
   """ This function is a generalisation of OPA2TPA and OPA23PA to arbitrary particle numbers.

      **PARAMETERS**
      OPAfreq:	frequencies of transitions in OPA. (Frequencies of modes * number of quanta in change)
	        array of lenght n
      freq00:	frequency of purely electronic transition
      OPAintens:intensities of the respective transitions in same order as OPAfreq
      intens00:	intensity of the purely electronic transition
	        array of lenght n
      mode:	number of vibrational states changing
	        array of lenght n

      **RETURNS**
      TPAfreq:	

      TPAfreq:	frequencies of the nPA-vibrational spectrum
      TPAintens:intensities of the nPA-vibrational spectrum	
   """
   def putN(j, n, intens, freq, mode, OPAintens, OPAfreq, oldmode):
      """ This function does the most calculation that is the iteration to the next number of particles
      """
      newintens=[]
      newfreq=[]
      #print 'arrived here:', mode, oldmode, n ,'-------------------'
      for i in range(len(intens)):
	 if intens[i]>0.000:
	    newintens.append(intens[i]) #this is OPA-part
	    newfreq.append(freq[i])
	    if n<=1:
	       continue 
	       #this saves creating new objects and running through loops without having results
	    tmpintens=[]
	    tmpfreq=[]
	    newmode=[]
	    nwemode=[]
	    for k in range(i+1, len(oldmode[0])): # go through whole range of states ...
	       tmpmode=[p for p in mode[:].T[i] if (p not in oldmode[:].T[k]) and np.all(p!=0) ]
	       #tmpmode=list(set(mode[:].T[i]).difference(oldmode[:].T[k]))
	       if tmpmode==[]:
		  continue
	       tempmode=oldmode[:].T[k]
	       if tempmode==[0]:
		  # that means, if mode[:].T[k] contains 0-0 transition
		  continue
	       try :
		  tmpmode=np.array(tmpmode[0])
		  #print tmpmode
		  xmode=[]
		  for j in range(len(tmpmode[0])):
		     xmode.append(tmpmode[0][j])
		  #### how to change rows and columns in a list?
		  newmode.append(xmode) #or xmode.T??
	       except IndexError:
		  newmode.append(float(tmpmode))
	       nwemode.append(tempmode[0])
	       #print "new modes:\n", newmode, nwemode, k
	       tmpintens.append(OPAintens[k]*intens[i])
	       tmpfreq.append(OPAfreq[k]+freq[i])
	    if len(tmpintens)>0:
	       xmode=[]
	       xmode.append(newmode)
	       if len(np.shape(xmode))>2:
		  xmode=np.matrix(xmode[0]).T
		  #print 'what I want:', xmode, np.shape(xmode)
		  nmode=np.zeros(( len(xmode)+1, len(xmode.T) ))
		  nmode[:-1]=xmode
		  nmode[-1]=nwemode
		  #print 'what I finally get:', nmode
	       else:
		  xmode.append(nwemode)
		  nmode=np.matrix(xmode)
	       #print "submitting:", nmode
	       #print tmpfreq
	       freq2, intens2=putN(i, n-1, tmpintens, tmpfreq, nmode, OPAintens, OPAfreq, oldmode)
	       for k in range(len(intens2)):
		  newintens.append(intens2[k])
		  newfreq.append(freq2[k])
      return np.array(newfreq), np.array(newintens)
	
   length=len(OPAfreq)
   for i in range(length):
      OPAfreq[i]-=freq00
      OPAintens[i]/=intens00
   newmode=np.zeros((1,len(mode))) #for n>1: matrix-structure needed
   newmode[0]=mode
   #print 'mode:', newmode
   TPAfreq, TPAintens=putN(-1, n, OPAintens, OPAfreq, newmode, OPAintens, OPAfreq, newmode)
   for i in range(len(TPAfreq)):
      TPAfreq[i]+=freq00
      TPAintens[i]*=intens00
   return TPAfreq, TPAintens

def OPA2TPA(OPAfreq,freq00, OPAintens,intens00, mode):
   length=len(OPAfreq)
   TPAfreq=np.zeros((length+1)*(length+2)//2+length+1) #this is overestimation of size...
   TPAintens=np.zeros((length+1)*(length+2)//2+length+1)
   TPAintens[0]=intens00 #this is OPA-part
   TPAfreq[0]=freq00
   ind=1
   for i in range(length):
      if mode[i]==0:
	 #0-0 transition is included, but only once! (before)
	 continue
      TPAintens[ind]=OPAintens[i] #this is OPA-part
      TPAfreq[ind]=OPAfreq[i]
      #print i, TPAfreq[ind]-freq00, TPAintens[ind]/intens00
      ind+=1
      for j in range(i+1,length):
	 ##not only same mode but all modes with lower number should not be taken into account here!?
	 if mode[i]==mode[j] or mode[j]==0: #both have same mode... or mode[j] is 0-0 transition
	    continue
	 TPAintens[ind]=OPAintens[i]*OPAintens[j]/intens00
	 TPAfreq[ind]=OPAfreq[i]+OPAfreq[j]-freq00
	 ind+=1
	 #print mode[i], mode[j]
   index=np.argsort(TPAfreq,kind='heapsort')
   TPAfreq=TPAfreq[index]
   TPAintens=TPAintens[index]
   #for i in range(len(TPAfreq)):
      #print TPAfreq[i], TPAintens[i]
   return TPAfreq, TPAintens

def OPA23PA(OPAfreq,freq00, OPAintens,intens00, mode):
   length=len(OPAfreq)
   TPAfreq=[]#np.zeros((((length+1)*(length+2)+3)*(length+3))//6+length+1) #this is overestimation of size...
   TPAintens=[]
   TPAintens.append(intens00) #this is OPA-part
   TPAfreq.append(freq00)
   for i in range(length):
      if mode[i]==0:
	 #0-0 transition is included, but only once!
	 continue
      TPAintens.append(OPAintens[i]) #this is OPA-part
      TPAfreq.append(OPAfreq[i])
      for j in range(i+1,length):
	 if mode[i]==mode[j] or mode[j]==0: #both have same mode... or mode[j] is 0-0 transition
	    continue
	 TPAintens.append(OPAintens[i]*OPAintens[j]/intens00)
	 TPAfreq.append(OPAfreq[i]+OPAfreq[j]-freq00)
	 for k in range(j+1,length):
	    if mode[k]==mode[j] or mode[j]==0: #both have same mode...
	       continue
	    if mode[k]==mode[i]:
	       continue
	    TPAintens.append(OPAintens[i]*OPAintens[j]*OPAintens[k]/(intens00*intens00))
	    TPAfreq.append(OPAfreq[i]+OPAfreq[k]+OPAfreq[j]-2*freq00)
   freq=np.zeros(len(TPAfreq))
   intens=np.zeros(len(TPAintens))
   for i in range(len(freq)): #this can not be done by np.matrix() due to dimensionality...
	freq[i]=TPAfreq[i]
	intens[i]=TPAintens[i]
   return freq, intens

def outspect(logging, T, opt, linspect, E=0):
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
   index=np.argsort(linspect[1], kind='heapsort')
   linspect[1]=linspect[1][index] #intensity
   linspect[2]=linspect[2][index] #mode
   linspect[0]=linspect[0][index] #frequency
   #find transition with minimum intensity to be respected
   minint=0
   for i in range(len(linspect[1])):
      if linspect[1][i]>=0.0001*linspect[1][-1]:
	 minint=i
	 break
   if logging[0]<3:
      logging[1].write('neglect '+repr(minint)+' transitions, use only '+repr(len(linspect[1])-minint)+" instead.\n")

      if logging[0]<2:
	 logging[1].write('minimal and maximal intensities:\n'+repr(linspect[1][minint])+' '+repr(linspect[1][-1])+"\n")

   #make nPA from OPA:
   if (re.search(r"to [\d]PA", opt, re.I) is not None) is True:
      n=re.findall(r"(?<=to )[\d](?=PA)", opt, re.I)
      if n[0]=='2':
	 ind=linspect[2].argmin()
	 #                          spectral frequency   0-0 transition      intensities          0-0 intensit.          modes        
	 TPAf, TPAi=OPA2TPA(linspect[0][minint:],linspect[0][ind] ,linspect[1][minint:], linspect[1][ind], 
											     linspect[2][minint:])
	 index=np.argsort(TPAi,kind='heapsort')
	 TPAi=TPAi[index] #resort by intensity
	 TPAf=TPAf[index]
	 minint=0
	 for i in range(len(TPAi)):
	    if TPAi[i]>=0.0001*TPAi[-1]:
	       minint=i
	       break
	 if logging[0]<3:
   	    logging[1].write('for TPA: again neglect '+repr(minint)+
		     ' transitions, use only '+repr(len(TPAi)-minint)+" instead.\n")
	 TPAintens=TPAi[minint:]
	 TPAfreq=TPAf[minint:]
      elif n[0]=='3':
	 ind=linspect[2].argmin()
	 TPAfreq, TPAintens=OPA2TPA(linspect[0][minint:],linspect[0][ind] ,linspect[1][minint:], 
			      linspect[1][ind], linspect[2][minint:])
	 minint=0
	 index=np.argsort(TPAintens,kind='heapsort')
	 TPAintens=TPAintens[index] #resort by intensity
	 TPAfreq=TPAfreq[index]
	 for i in range(len(TPAintens)):
	    if TPAintens[i]>=0.0001*TPAintens[-1]:
	       minint=i
	       break
	 if logging[0]<3:
	    logging[1].write('for 3PA: again neglect '+repr(minint)+
		     ' transitions, use only '+repr(len(TPAintens)-minint)+" instead.\n")
	 TPAintens=TPAintens[minint:] #resort by intensity
	 TPAfreq=TPAfreq[minint:]
      elif n[0]==1:
	 TPAfreq=linspect[0][minint:]
	 TPAintens=linspect[1][minint:]
      else:
	 logging[1].write("to <n>PA was given but not recognised.\n")
   else:
      n=re.findall(r"(?<=to nPA:)[ \d]*", opt, re.I)
      if n==[]:
	 TPAfreq=linspect[0][minint:]
	 TPAintens=linspect[1][minint:]
      else:
	 n=float(n[0])
	 ind=linspect[2].argmin()
	 #                          spectral frequency   0-0 transition      intensities          0-0 intensit.          modes        
	 TPAfreq, TPAintens=OPA2nPA(linspect[0][minint:],linspect[0][ind] ,linspect[1][minint:], linspect[1][ind], 
											     linspect[2][minint:], n)
	 index=np.argsort(TPAintens,kind='heapsort')
	 TPAintens=TPAintens[index] #resort by intensity
	 TPAfreq=TPAfreq[index]
	 minint=0
	 for i in range(len(TPAintens)):
	    if TPAintens[i]>=0.00005*TPAintens[-1]:
	       minint=i
	       break
	 if logging[0]<3:
	    logging[1].write('for {0}PA: again neglect {1} \n'.format(n, minint)+
		     ' transitions, use only '+repr(len(TPAintens)-minint)+" instead.\n")
	 TPAintens=TPAintens[minint:] #resort by intensity
	 TPAfreq=TPAfreq[minint:]


   #find transition with minimum intensity to be respected
   #the range of frequency ( should be greater than the transition-frequencies)
   if omega==None:
      if minfreq==0:
	 minfreq=np.min(TPAfreq)-20-gamma
      if maxfreq==0:
	 maxfreq=np.max(TPAfreq)+20+gamma
   else:
      minfreq=omega[0]
      maxfreq=omega[-1]
   if logging[0]<3:
      logging[1].write('maximal and minimal frequencies:\n {0} {1}\n'.format(maxfreq, minfreq))
   #truncate arrays and sort by index for further efficient processes
   #if no other grid is defined: use linspace in range
   if omega==None:
      omega=np.linspace(minfreq,maxfreq,gridpt)
      if logging[0]<2:
	 logging[1].write("omega is equally spaced\n")
   spect=np.zeros(len(omega))
   sigma=gamma*2/2.355 #if gaussian used: same FWHM

   index=np.argsort(TPAfreq,kind='heapsort') #sort by freq
   freq=TPAfreq[index]
   intens=TPAintens[index]
   if logging[0]<1:
      logging[1].write('intensity   frequency\n')
      for i in range(len(intens)):
	 logging[1].write(u"{0} {1}\n".format(intens[i], freq[i]))

   mini=0
   maxi=len(freq) #just in case Gamma is too big or frequency-range too low
   for i in range(1,len(freq)):
      if freq[i]>=5*gamma+freq[0]:
	 maxi=i
	 break
   if spectfile==None:
      spectfile="/dev/null" #discart spectrum
   out = open(spectfile, "w")
   logging[1].write("broadened spectrum:\n frequency      intensity\n")
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
	 #out.write(u" {0}  {1}\n".format(omega[i] ,spect[i]))
	 logging[1].write(u" {0}  {1}\n".format(omega[i] ,spect[i]))
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
	 logging[1].write(u" {0}  {1}\n".format(omega[i] ,spect[i]))
   out.close()

version=1.1
# End of broadening.py
