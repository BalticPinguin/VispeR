#!/usr/bin/python
# filename: broadening.pyx
import numpy as np, re

# Below are the conversion factors and fundamental constant

def handel_input(opt):
   #set default values (to have all variables set)
   gridfile=None
   gamma=1 #by default: only slight broadening
   gridpt=5000
   omega=None
   minfreq=0
   maxfreq=0
   shape='g'
   stick=False

   tmpgrid=re.findall(r"(?<=grid=)[ \=\s\w\.;]+", opt, re.M)
   if len(tmpgrid)==1: 
      # i.e. if grid is specified
      grid=re.findall(r"[\w\.]+", tmpgrid[0], re.M)
      if len(grid)==1:
         #that means, if either one number (# of gridpoints or a file) is given
         try:
            gridpt=float(grid[0])
         except ValueError: # if grid is no a number
            gridfile=grid[0]
      elif len(grid)==3:
         # that means there is the number of gridpoints, min- and max frequency given
         gridpt=float(grid[0])
         minfreq=float(grid[1])
         maxfreq=float(grid[2])
      if gridfile!=None:
         #read file in format of linspect
         grid=[]
         with open(gridfile) as f:
            lis=[line.split() for line in f]  # create a list of lists
            for i,x in enumerate(lis):        # get the list items 
               grid.append(float(x[0]))
         omega=np.zeros(len(grid))
         for i in range(len(grid)):
            omega[i]=grid[i]
   #see, whether a broadening is given
   if (re.search(r"(?<=gamma=)[ \d\.,]+", opt, re.I) is not None)  is True:
      gamma=re.findall(r"(?<=gamma=)[ \d\.]+", opt, re.I)
      gamma=float(gamma[0])

   shape=re.findall(r"(?<=shape=)[ \w]+", opt, re.I)
   if len(shape)>0:
      # there are several options each
      if shape[0] in ["lorentzian", "Lorentzian", "L", "l"]:
         shape="l"
      elif shape[0] in ["gaussian", "Gaussian", "G", "g"]:
         shape="g"

   if (re.search(r"stick", opt, re.I) is not None) is True:
      stick=True

   spectfile=re.findall(r"(?<=spectfile=)[\w.]+", opt, re.I)
   if spectfile==[]:
      spectfile=re.findall(r"(?<=spectfile= )[\w.]+", opt, re.I)
      if spectfile==[]:
         spectfile=None
      else:
         spectfile=spectfile[-1]
   else:
      spectfile=spectfile[-1]
   return omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape, stick

def OPA2nPA(logwrite, OPAfreq,freq00, OPAintens, intens00, mode, n, stick):
   """ This function is a generalisation of OPA2TPA and OPA23PA to arbitrary particle numbers.

   **PARAMETERS**
   logging: This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
            and 4- only main information) and logging[1] is the file, already opened, to write the information in.
   OPAfreq:  frequencies of transitions in OPA. (Frequencies of modes * number of quanta in change)
             array of lenght n
   freq00:   frequency of purely electronic transition
   OPAintens:intensities of the respective transitions in same order as OPAfreq
   intens00: intensity of the purely electronic transition
             array of lenght n
   mode:     number of vibrational states changing
             array of lenght n
   n:        number of maximal changing vibrations simultanously
   stick:    a boolean variable, stating whether to print the stick-spectrum or not

   **RETURNS**
   TPAfreq:  

   TPAfreq:  frequencies of the nPA-vibrational spectrum
   TPAintens:intensities of the nPA-vibrational spectrum     
   """
   append=np.append

   def allnonzero(foo):
      """ a more efficient version of np.all(foo!=0), made for arrays and integers...
      """
      cdef int s
      try:
         foo=np.array(foo)
         for s in range(len(foo)):
            if foo[s]==0:
               return False
      except TypeError:
         if foo==0:
            return False
      return True
     
   def putN(int j, int n, intens, freq, mode, OPAintens, OPAfreq, oldmode, stick, logwrite):
      """ This function does the most calculation that is the iteration to the next number of particles
          therefor it is highly optimised into C
      """
      cdef int i
      cdef int k
      cdef double intensi
      cdef double freqi
      newintens=[]
      newfreq=[]
      
      for i in xrange(len(intens)):
         intensi=intens[i]
         if intensi>1e-9: # only use this if the intensity is reasonably high
            newintens=append(newintens, intensi) #this is OPA-part
            freqi=freq[i]
            newfreq=append(newfreq, freqi)
            if stick:
               logwrite(u" %s  %s\n" %(intensi, freqi))

            if n<=1:
               continue 
               #this saves creating new objects and running through loops without having results
            tmpintens=[]
            tmpfreq=[]
            newmode=[]
            nwemode=[]
            for k in range(i+1, len(oldmode[0])): # go through whole range of states ...
               tempmode=oldmode[:].T[k]
               if tempmode==[0]:
                  # that means, if mode[:].T[k] contains 0-0 transition
                  continue
               #tmpmode=np.array(mode[:].T[i])
               tmpmode=mode[:].T[i]
               if tmpmode==[]:
                  continue
               if not allnonzero(tmpmode):
                  continue
               if tempmode in tmpmode:
                  continue
               foo=newmode
               foo.append(tmpmode)
               newmode=foo
               nwemode=append(nwemode, tempmode)
               tmpintens=append(tmpintens, OPAintens[k]*intensi)
               tmpfreq=append(tmpfreq, OPAfreq[k]+freqi)
            if len(tmpintens)>0:
               xmode=newmode
               if np.shape(xmode)[1]>=2:
                  xmode=np.matrix(xmode).T
                  nmode=np.zeros(( len(xmode)+1, len(xmode.T) ))
                  nmode[:-1]=xmode
                  nmode[-1]=nwemode
               else:
                  nmode=np.zeros(( 2 , len(xmode) ))
                  nmode[0]=xmode
                  nmode[1]=nwemode
               freq2, intens2=putN(i, n-1, tmpintens, tmpfreq, nmode, OPAintens, OPAfreq, oldmode, stick, logwrite)
               for k in range(len(intens2)):
                  newintens=append(newintens, intens2[k])
                  newfreq=append(newfreq, freq2[k])
      return np.array(newfreq), np.array(newintens)
        
   length=len(OPAfreq)
   for i in range(length):
      OPAfreq[i]-=freq00
      OPAintens[i]/=intens00
   newmode=np.zeros((1,len(mode))) #for n>1: matrix-structure needed
   newmode[0]=mode
   x=mode.max()
   if n>x:
      n=x
      #save complexity, that it does not try to make more combinations than actually possible...
   #np.set_printoptions(precision=5, linewidth=138)
   if stick:
      logwrite("stick-spectrum in n-particle approximation\n")
      logwrite(" intensity    frequency \n")
   TPAfreq, TPAintens=putN(-1, n, OPAintens, OPAfreq, newmode, OPAintens, OPAfreq, newmode, stick, logwrite)
   for i in xrange(len(TPAfreq)):
      TPAfreq[i]+=freq00
      TPAintens[i]*=intens00
   return TPAfreq, TPAintens

def OPA2TPA(logwrite, OPAfreq,freq00, OPAintens,intens00, mode, stick):
   """ This function calculates  a Two-particle spectra using one-particle spectra.

   **PARAMETERS**
   logging: This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
            and 4- only main information) and logging[1] is the file, already opened, to write the information in.
   OPAfreq:  frequencies of transitions in OPA. (Frequencies of modes * number of quanta in change)
             array of lenght n
   freq00:   frequency of purely electronic transition
   OPAintens:intensities of the respective transitions in same order as OPAfreq
   intens00: intensity of the purely electronic transition
             array of lenght n
   mode:     number of vibrational states changing
             array of lenght n
   stick:    a boolean variable, stating whether to print the stick-spectrum or not

   **RETURNS**
   TPAfreq:  

   TPAfreq:  frequencies of the nPA-vibrational spectrum
   TPAintens:intensities of the nPA-vibrational spectrum     
   """
   length=len(OPAfreq)
   TPAfreq=np.zeros((length+1)*(length+2)//2+length+1) #this is overestimation of size...
   TPAintens=np.zeros((length+1)*(length+2)//2+length+1)
   if stick:
      logwrite(u"stick-spectrum in 3-particle approximation:\n  intensity   frequency")
   ind=0
   for i in range(length):
      TPAintens[ind]=OPAintens[i] #this is OPA-part
      TPAfreq[ind]=OPAfreq[i]
      ind+=1
      if stick:
         logwrite(u" %f  %e\n"%(OPAintens[i], OPAfreq[i]))
      for j in range(i+1,length):
         ##not only same mode but all modes with lower number should not be taken into account here!?
         if mode[i]==mode[j] or mode[j]==0: #both have same mode... or mode[j] is 0-0 transition
            continue
         if OPAintens[i]*OPAintens[j]/intens00<intens00*0.0001:
            #save memory by not saving low-intensity-modes
            continue
         TPAintens[ind]=OPAintens[i]*OPAintens[j]/intens00
         TPAfreq[ind]=OPAfreq[i]+OPAfreq[j]-freq00
         ind+=1
         if stick:
            logwrite(u" %f  %e\n"%(TPAintens[-1], TPAfreq[-1]))
   index=np.argsort(TPAfreq,kind='heapsort')
   TPAfreq=TPAfreq[index]
   TPAintens=TPAintens[index]
   return TPAfreq, TPAintens

def OPA23PA(logwrite, OPAfreq,freq00, OPAintens,intens00, mode, stick):
   """ This function calculates  a Three-particle spectra using one-particle spectra.

   **PARAMETERS**
   logging: This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
            and 4- only main information) and logging[1] is the file, already opened, to write the information in.
   OPAfreq:  frequencies of transitions in OPA. (Frequencies of modes * number of quanta in change)
             array of lenght n
   freq00:   frequency of purely electronic transition
   OPAintens:intensities of the respective transitions in same order as OPAfreq
   intens00: intensity of the purely electronic transition
             array of lenght n
   mode:     number of vibrational states changing
             array of lenght n
   stick:    a boolean variable, stating whether to print the stick-spectrum or not

   **RETURNS**
   TPAfreq:  

   TPAfreq:  frequencies of the nPA-vibrational spectrum
   TPAintens:intensities of the nPA-vibrational spectrum     
   """
   length=len(OPAfreq)
   TPAfreq=[]#np.zeros((((length+1)*(length+2)+3)*(length+3))//6+length+1) #this is overestimation of size...
   TPAintens=[]
   if stick:
      logwrite(u"stick-spectrum in 3-particle approximation:\n  intensity   frequency")
   TPAintens.append(intens00) #this is OPA-part
   TPAfreq.append(freq00)
   # go through the whole spectrum (besides 0-0) and compute all combinations besides self-combination
   if stick:
      logwrite(u" %f  %e\n"%(TPAintens[0], TPAfreq[0]))
   for i in range(length):
      if mode[i]==0:
         #0-0 transition is included, but only once!
         continue
      TPAintens.append(OPAintens[i]) #this is OPA-part
      TPAfreq.append(OPAfreq[i])
      if stick:
         logwrite(u" %f  %e\n"%(OPAintens[i], OPAfreq[i]))
      # here the combination part starts
      for j in range(i+1,length):
         if mode[i]==mode[j] or mode[j]==0: #both have same mode... or mode[j] is 0-0 transition
            continue
         if OPAintens[i]*OPAintens[j]/intens00<intens00*0.0002:
            #save memory by not saving low-intensity-modes
            continue
         TPAintens.append(OPAintens[i]*OPAintens[j]/intens00)
         TPAfreq.append(OPAfreq[i]+OPAfreq[j]-freq00)
         if stick:
            logwrite(u" %f  %e\n"%(TPAintens[-1], TPAfreq[-1]))
         # look for all combinations of combinations for three-particle approx.
         for k in range(j+1,length):
            if mode[k]==mode[j] or mode[k]==0: #both have same mode...
               continue
            if mode[k]==mode[i]:
               continue
            if OPAintens[i]*OPAintens[j]*OPAintens[k]/(intens00*intens00)<intens00*0.0001:
               #save memory by not saving low-intensity-modes
               continue
            TPAintens.append(OPAintens[i]*OPAintens[j]*OPAintens[k]/(intens00*intens00))
            TPAfreq.append(OPAfreq[i]+OPAfreq[k]+OPAfreq[j]-2*freq00)
            if stick:
               logwrite(u" %f  %e\n"%(TPAintens[-1], TPAfreq[-1]))
   # save the spectrum into numpy-matrices
   freq=np.zeros(len(TPAfreq))
   intens=np.zeros(len(TPAintens))
   for i in xrange(len(freq)): #this can not be done by np.matrix() due to dimensionality...
        freq[i]=TPAfreq[i]
        intens[i]=TPAintens[i]
   return freq, intens

def concise(intens, freq, broadness):
   """
   This function shrinks length of the stick-spectrum to speed-up the calculation of broadened spectrum 
   (folding with lineshape-function).
   It puts all transitions within a tenth of the Gaussian-width into one line.

   ==PARAMETERS==
   intens:      intensity (of stick-spectrum-points) sorted by increasing freq
   freq:        frequency (of stick-spectrum-points) sorted by increasing freq
   broadness:   gamma from the Lorentian-courve; specifying, how many lines will be put together

   ==RETURNS==
   intens2:     shrinked intensity-vector, sorted by increasing frequency
   freq2:       shrinked frequency-vector, sorted by increasing frequency
   """
   # both arrays are frequency-sorted 
   #initialize arrays
   intens2=[]
   freq2=[]
   mini=0
   #go through spectrum and sum everything up.
   for i in range(len(freq)):
      # index-range, put into one line is broadness/5; this should be enough
      tempintens=0
      for j in xrange(mini, len(freq)):
         tempintens+=intens[j]
         if freq[j]>=freq[mini]+broadness/5.:
            mini=j # set mini to its new value
            intens2.append(tempintens)
            freq2.append(freq[j]-broadness/10.)
            break
   return intens2, freq2

def outspect(logging, float T, opt, linspect, float E=0):
   """This function calculates the broadened spectrum given the line spectrum, 
   frequency-rage and output-file whose name is first argument. 
   As basis-function a Lorentzian is assumed with a common width.
   
   **PARAMETERS:**
   logging: This variable consists of two parts: logging[0] specifies the level of print-out (which is between 0- very detailed
            and 4- only main information) and logging[1] is the file, already opened, to write the information in.
   T:       temperature of the system
   opt:     a string that contains all options that were given for this part in the input-file. See documentation 
            for details of it's allowed/used content
   linspect:The line-spectrum that has to be broadened: A array/matrix with 3(!!) rows: 
            Frequency, intentensity and mode number (last one is important for making multiple-particle spectra 
   E:       energy-shift of the 0-0 transition. Important if the excited state is not the lowest and
            thermal equilibration with the lower states should be considered

   **RETURNS:**
   nothing; the key values (broadened spectra/ many-particle-app. linespectra) are printed into log-files.
   
   """
   cdef float minint=0
   cdef int i

   omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape , stick=handel_input(opt)
   #read file in format of linspect

   #sort spectrum with respect to size of elements
   index=np.argsort(linspect[1], kind='heapsort')
   linspect[1]=linspect[1][index] #intensity
   linspect[2]=linspect[2][index] #mode
   linspect[0]=linspect[0][index] #frequency
   #find transition with minimum intensity to be respected

   #truncate all transitions having less than 0.0001% of
   for i in range(len(linspect[1])):
      if linspect[1][i]>=1e-6*linspect[1][-1]:
         minint=i
         break
   if logging[0]<3:
      logging[1].write('neglect '+repr(minint)+' transitions, use only '+repr(len(linspect[1])-minint)+" instead.\n")

      if logging[0]<2:
         logging[1].write('minimal and maximal intensities:\n'+repr(linspect[1][minint])+' '+repr(linspect[1][-1])+"\n")

   logwrite=logging[1].write  #important for later loops: avoiding '.'s speeds python-codes up!!
   #make nPA from OPA:
   if (re.search(r"to [\d]PA", opt, re.I) is not None) is True:
      n=re.findall(r"(?<=to )[\d](?=PA)", opt, re.I)
      if n[0]=='2':
         ind=linspect[2].argmin()
         #  spectral frequency   0-0 transition   intensities      0-0 intensit.          modes        
         TPAf, TPAi=OPA2TPA(logwrite, linspect[0][minint:],linspect[0][ind] ,
                            linspect[1][minint:], linspect[1][ind], linspect[2][minint:], stick)
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
         TPAfreq, TPAintens=OPA23PA(logwrite, linspect[0][minint:],linspect[0][ind] ,linspect[1][minint:], 
                              linspect[1][ind], linspect[2][minint:], stick)
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
         if stick:
            logwrite=logging[1].write
            logwrite("stick-spectrum in one-particle approximation")
            logwrite(" intensity    frequency ")
            for k in range(len(TPAfreq)):
               logwrite(u" %s  %s\n" %(TPAintens[k], TPAfreq[k]))
      else:
         TPAfreq=linspect[0][minint:]
         TPAintens=linspect[1][minint:]
         logging[1].write("to <n>PA was given but not recognised.\n")
         if stick:
            logwrite=logging[1].write
            logwrite("stick-spectrum in one-particle approximation")
            logwrite(" intensity    frequency ")
            for k in range(len(TPAfreq)):
               logwrite(u" %s  %s\n" %(TPAintens[k], TPAfreq[k]))
   else:
      n=re.findall(r"(?<=to nPA:)[ \d]*", opt, re.I)
      if n==[]:
         TPAfreq=linspect[0][minint:]
         TPAintens=linspect[1][minint:]
         if stick:
            logwrite=logging[1].write
            logwrite("stick-spectrum in one-particle approximation")
            logwrite(" intensity    frequency ")
            for k in range(len(TPAfreq)):
               logwrite(u" %s  %s\n" %TPAintens[k] %TPAfreq[k])
      else:
         n=float(n[0])
         ind=linspect[2].argmin() #get position of 0-0 transition
         #                          spectral frequency   0-0 transition      intensities          0-0 intensit.          modes        
         TPAfreq, TPAintens=OPA2nPA(logwrite, linspect[0][minint:],linspect[0][ind] ,
                     linspect[1][minint:], linspect[1][ind], linspect[2][minint:], n, stick)
         index=np.argsort(TPAintens,kind='heapsort')
         TPAintens=TPAintens[index] #resort by intensity
         TPAfreq=TPAfreq[index]
         minint=0
         for i in range(len(TPAintens)):
            if TPAintens[i]>=1e-6*TPAintens[-1]:
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
         minfreq=np.min(TPAfreq)-20-gamma*15
      if maxfreq==0:
         maxfreq=np.max(TPAfreq)+20+gamma*15
   else:
      minfreq=omega[0]
      maxfreq=omega[-1]
   if logging[0]<3:
      logging[1].write('maximal and minimal frequencies:\n {0} {1}\n'.format(maxfreq, minfreq))
   #if no other grid is defined: use linspace in range
   if omega==None:
      omega=np.linspace(minfreq,maxfreq,gridpt)
      if logging[0]<2:
         logging[1].write("omega is equally spaced\n")
   sigma=gamma*2/2.355 #if gaussian used: same FWHM

   #TPAintens/=TPAintens[-1]*0.01 #normalise spectrum (due to highest peak) to make different approaches comparable

   index=np.argsort(TPAfreq,kind='heapsort') #sort by freq
   freq=TPAfreq[index]
   intens=TPAintens[index]


   if logging[0]<1:
      logwrite('truncated and sorted stick-spectrum '
                     '(this is what will be taken into account for broadening):\n')
      logwrite('intensity   frequency\n')
      for i in range(len(intens)):
         logwrite(u"%f   %f\n"%(intens[i], freq[i]))
   mini=0
   if spectfile==None:
      out=logging[1]
   else:
      out = open(spectfile, "w")

   spect=np.zeros(len(omega))
   if spectfile==None: #that means spectrum is printed into log-file
      logwrite("broadened spectrum:\n frequency      intensity\n")
   outwrite=out.write
   #this shrinks the size of the spectral lines; hopefully accelerates the script.
   #intens, freq=concise(intens,freq, sigma)

   maxi=len(freq)-1 #just in case Gamma is too big or frequency-range too low
   for i in range(0,len(freq)-1):
      if freq[i]>=10*gamma+freq[0]:
         maxi=i
         break
   if shape=='g':
      sigmasigma=2.*sigma*sigma
      npexp=np.exp
      for i in xrange(len(omega)): 
         omegai=omega[i]
         for j in range(maxi,len(freq)-1):
            if freq[j]>=10*gamma+omegai:
               maxi=j
               break
         for j in range(max(0,mini),maxi):
            if freq[j]>=omegai-8*gamma:
               mini=max(j-1,0) # else it becomes -1 and hence the spectrum is wrong
               break
         for k in range(mini,min(maxi,len(freq))):
            spect[i]+=intens[k]*npexp(-(omegai-freq[k])*(omegai-freq[k])/(sigmasigma))
         outwrite(u" %f  %e\n" %(omegai, spect[i]))
   else:  #shape=='l':
      gammagamma=gamma*gamma
      for i in xrange(len(omega)): 
         omegai=omega[i]
         for j in xrange(maxi,len(freq)-1):
            if freq[j]>=20*gamma+omegai:
               maxi=j
               break
         else: #if it never reached 'break'-statement
            maxi=len(freq) 
         for j in xrange(max(0,mini),maxi):
            if freq[j]>=omegai-20*gamma:
               mini=max(j-1,0)
               break
               #this should always be reached somewhen
         for k in range(mini,min(maxi,len(freq))):
            spect[i]+=intens[k]/((omegai-freq[k])*(omegai-freq[k])+gammagamma)
         outwrite(u" %f   %e\n" %(omegai, spect[i]))
   if spectfile!=None:
      #only close file if it was opened here
      out.close()

version=1.5
# End of broadening.py
