#!/usr/bin/python2
# filename: output_spect.py

#include [[MultiPart.py]]
import MultiPart
import numpy as np
import re 

# CHANGELOG
# =========
#in version 0.1.0:  
#   1) intitialised class
#   2) fixed class-related issues (missing self.)
#   3) Added function concise to speed-up calculations.
#

class broaden():
   """ This class computes the broadened spectrum out of a stick-spectrum.
         
       **Interface functions **
       init     - initialisation of the class. Its argument is a pointer 
                  to a class of type Spect (or inherited from it); than
                  all required quantities are copied from there.
       outspect - computes the broadened spectrum and saves it to a file.
   """
   #set default values (to have all variables set)
   gridfile=None
   gamma=1 #by default: only slight broadening
   gridpt=5000
   omega=None
   minfreq=0
   maxfreq=0
   shape='g'
   stick=False
   log=""

   #BEGIN DEFINITION OF METHODS
   def __init__(self, parent):
      """ initialise the quantities needed for the calculation of
         broadened spectrum
      """
      #get a pointer to parent-object that e.g. has the stick-spectrum
      self.parent=parent
      self.log=self.parent.log

      #get some information about the grid
      tmpgrid=re.findall(r"(?<=grid=)[ \=\s\w\.;]+", parent.broadopt, re.M)
      if len(tmpgrid)==1: 
         # i.e. if grid is specified
         grid=re.findall(r"[\w\.]+", tmpgrid[0], re.M)
         if len(grid)==1:
            #that means, if either one number (# of gridpoints or a file) is given
            try:
               self.gridpt=float(grid[0])
            except ValueError: # if grid is no a number
               self.gridfile=grid[0]
         elif len(grid)==3:
            # that means there is the number of gridpoints, min- and max frequency given
            self.gridpt=float(grid[0])
            self.minfreq=float(grid[1])
            self.maxfreq=float(grid[2])
   
         if self.gridfile!=None:
               #read file in format of spect
               grid=[]
               with open(self.gridfile) as f:
                  lis=[line.split() for line in f]  # create a list of lists
                  for i,x in enumerate(lis):        # get the list items 
                     grid.append(float(x[0]))
               self.omega=np.zeros(len(grid))
               for i in range(len(grid)):
                  self.omega[i]=grid[i]

      #see, whether a broadening is given
      if (re.search(r"(?<=gamma=)[ \d\.,]+", parent.broadopt, re.I) is not None)  is True:
         gamma=re.findall(r"(?<=gamma=)[ \d\.]+", parent.broadopt, re.I)
         self.gamma=float(gamma[0])
   
      #which line shape function should be used?
      # FIXME: also some Voigt-profile should be possible here.
      shape=re.findall(r"(?<=shape=)[ \w]+", parent.broadopt, re.I)
      if len(shape)>0:
      # there are several options each
         if shape[0] in ["lorentzian", "Lorentzian", "L", "l"]:
            self.shape="l"
         elif shape[0] in ["gaussian", "Gaussian", "G", "g"]:
            self.shape="g"

      # should the stick-spectrum be printed?
      if (re.search(r"stick", parent.broadopt, re.I) is not None) is True:
         self.stick=True
     
      #output-file for the broadened spectrum
      spectfile=re.findall(r"(?<=spectfile=)[\w._\-]+", parent.broadopt, re.I)
      if spectfile==[]:
         spectfile=re.findall(r"(?<=spectfile= )[\w._\-]+", parent.broadopt, re.I)
         if spectfile==[]:
            self.spectfile=None
         else:
            self.spectfile=spectfile[-1]
      else:
         self.spectfile=spectfile[-1]

   def outspect(self, E=0):
      """This function calculates the broadened spectrum given the line spectrum, 
         frequency-rage and output-file whose name is first argument. 
         As basis-function a Lorentzian is assumed with a common width.
      
         **PARAMETERS:**
         E:       energy-shift of the 0-0 transition. Important if the excited 
                  state is not the lowest and
                  thermal equilibration with the lower states should be considered
   
         **RETURNS:**
         nothing; the key values (broadened spectra/ many-particle-app. 
                  linespectra) are printed into log-files.
         
      """

      minint=0
      self.log.write("\n STARTING TO CALCULATE BROADENED SPECTRUM.\n")
      self.spect=self.parent.spect

      #sort spectrum with respect to size of elements
      index=np.argsort(self.spect[1], kind='heapsort')
      self.spect[0]=self.spect[0][index] #frequency
      self.spect[1]=self.spect[1][index] #intensity
      self.spect[2]=self.spect[2][index] #mode
      #find transition with minimum intensity to be respected
   
      #truncate all transitions having less than 0.0001% of
      for i in range(len(self.spect[1])):
         if self.spect[1][i]>=1e-6*self.spect[1][-1]:
            minint=i
            break
      self.log.write('neglect '+repr(minint)+' transitions, use only '+
                                repr(len(self.spect[1])-minint)+" instead.\n", 3)
   
      self.log.write('minimal and maximal intensities:\n'+
              repr(self.spect[1][minint])+' '+repr(self.spect[1][-1])+"\n", 2)
      
      #important for later loops: avoiding '.'s speeds python-codes up!!
      logwrite=self.log.write  

      #make nPA from OPA if requested.
      n=re.findall(r"(?<=to nPA:)[ \d]*", self.parent.broadopt, re.I)
      if n!=[]:
         MakeFull=MultiPart.OPAtoNPA(float(n[-1].strip()))
         self.log.write("\n REACHING OPA TO NPA-PART. \n")
         self.log.write(" ----------------------------------------"+
                                                 "-------- \n")
         MakeFull.GetSpect(self.spect, minint)
         TPAintens, TPAfreq=MakeFull.Calc()
      else: 
         TPAfreq=self.spect[0][minint:]
         TPAintens=self.spect[1][minint:]

      #print the stick-spectrum to a .stick-file, if it was requested:
      if self.stick:
         stickfile=self.parent.log.logfile.split(".")[0]
         stickout=open(stickfile+".stick", "w")
         stickout.write(" Intensity  \t frequency \n")
         for i in range(len(TPAfreq)):
            stickout.write(" %3.6g  \t %3.6f\n"%(TPAintens[i],TPAfreq[i]))
         stickout.close()
            
      #find transition with minimum intensity to be respected
      #the range of frequency ( should be greater than the transition-frequencies)
      if self.omega==None:
         if self.minfreq==0:
            self.minfreq=np.min(TPAfreq)-20-self.gamma*15
         if self.maxfreq==0:
            self.maxfreq=np.max(TPAfreq)+20+self.gamma*15
      else:
         self.minfreq=self.omega[0]
         self.maxfreq=self.omega[-1]
      self.log.write('maximal and minimal frequencies:\n {0} {1}\n'.format(self.maxfreq, self.minfreq), 3)
      #if no other grid is defined: use linspace in range
      if self.omega==None:
         self.omega=np.linspace(self.minfreq,self.maxfreq,self.gridpt)
         self.log.write("omega is equally spaced\n",2)
   
      
      if self.gamma*1.1<=(self.maxfreq-self.minfreq)/self.gridpt:
         self.log.write("\n WARNING: THE GRID SPACING IS LARGE COMPARED TO THE WIDTH OF THE PEAKS.\n"
              "THIS CAN ALTER THE RATIO BETWEEN PEAKS IN THE BROADENED SPECTRUM!")
   
      index=np.argsort(TPAfreq,kind='heapsort') #sort by freq
      freq=TPAfreq[index]
      intens=TPAintens[index]
   
      if self.spectfile==None:
         self.out=self.log
      else:
         self.out = open(self.spectfile, "w")
   
      if self.spectfile==None: #that means spectrum is printed into log-file
         logwrite("broadened spectrum:\n frequency      intensity\n")
      if self.shape=='g':
         self.__broadenGauss(intens, freq)
      else:  #shape=='l':
         self.__broadenLorentz(intens, freq)
      if self.spectfile!=None:
         #only close file if it was opened here
         self.out.close()

   def __broadenGauss(self, intens, freq):
      """This private function does the actual computation of the spectrum
         in case of a Gaussian as line spape function.
         The intensities and frequencies are required to be sorted by increasing
         frequencies.
      """
      sigma=self.gamma*2/2.355 #if gaussian used: same FWHM
      omega=self.omega
      sigmasigma=2.*sigma*sigma # these two lines are to avoid multiple calculations of the same
      npexp=np.exp
      intens/=sigma # scale all transitions according to their width.
      outwrite=self.out.write

      #this shrinks the size of the spectral lines; hopefully accelerates the script.
      intens, freq=concise(intens, freq, sigma)
      lenfreq=len(freq)
      maxi=lenfreq-1 #just in case Gamma is too big or frequency-range too low
      mini=0
      # set the largest index to be taken into account for the first transition.
      for i in range(0,lenfreq-1):
         if freq[i]>=10*sigma+freq[0]:
            maxi=i
            break

      # go through grid for the broadened spectrum and
      # compute the intensity at this point.
      for i in xrange(len(omega)): 
         omegai=omega[i]
         #re-adjust the limits for the transitions to be taken into account at this point
         for j in range(maxi,lenfreq):
            if freq[j]>=10*sigma+omegai:
               maxi=j
               break
         for j in range(mini,maxi):
            if freq[j]>=omegai-10*sigma:
               # else it becomes -1 and hence the spectrum is wrong
               mini=max(j-1,0) 
               break
         spect=0
         #sum up all contributions of neighbouring transitions.
         for k in range(mini,maxi+1):
            spect+=intens[k]*npexp(-(omegai-freq[k])*(omegai-freq[k])/(sigmasigma))
         #write the value to file
         outwrite(u" %f  %e\n" %(omegai, spect))

   def __broadenLorentz(self, intens, freq):
      """This private function does the actual computation of the spectrum
         in case of a Lorentzian as line spape function.
         The intensities and frequencies are required to be sorted by increasing
         frequencies.
      """
      gamma=self.gamma
      omega=self.omega
      gammagamma=gamma*gamma
      outwrite=self.out.write
      
      #this shrinks the size of the spectral lines; hopefully accelerates the script.
      intens, freq=concise(intens,freq, sigma)
      lenfreq=len(freq)
      maxi=lenfreq-1 #just in case Gamma is too big or frequency-range too low
      mini=0
      for i in range(0,lenfreq-1):
         if freq[i]>=10*self.gamma+freq[0]:
            maxi=i
            break
      
      # go through grid for the broadened spectrum and
      # compute the intensity at this point.
      for i in xrange(len(omega)): 
         omegai=omega[i]
         #re-adjust the limits for the transitions to be taken into account at this point
         for j in range(maxi,lenfreq):
            if freq[j]>=30*gamma+omegai:
               maxi=j
               break
         for j in range(mini,maxi):
            if freq[j]>=omegai-30*gamma:
               # else it becomes -1 and hence the spectrum is wrong
               mini=max(j-1,0) 
               break
         omegai=omega[i]
         spect=0
         #sum up all contributions of neighbouring transitions.
         for k in range(mini,maxi+1):
            spect+=intens[k]*gamma/((omegai-freq[k])*(omegai-freq[k])+gammagamma)
         #write the value to file
         outwrite(u" %f   %e\n" %(omegai, spect))
   #END DEFINITION OF METHODS
   
def concise(intens, freq, sigma):
   """ This function shrinks length of the stick-spectrum to speed-up the 
      calculation of broadened spectrum (folding with lineshape-function).
      It puts all transitions within a tenth of the Gaussian-width into one line.
   
      ==PARAMETERS==
      broadness:   gamma from the Lorentian-courve; specifying, 
                  how many lines will be put together
   
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
         if freq[j]>=freq[mini]+sigma/5.:
            mini=j # set mini to its new value
            intens2.append(tempintens)
            freq2.append(freq[j]-sigma/10.)
            break
   return intens2, freq2

version='0.1'
#End of output_spect.py
