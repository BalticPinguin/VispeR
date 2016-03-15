#!/usr/bin/python2
# filename: Spect.py
import numpy as np
from scipy.linalg.lapack import dsyev
import re, mmap, os.path, math, sys
import Read
import MultiPart
import random

# CHANGELOG
# =========
#in version 0.1.6:  
#   a) added RMSD_reorient   
#   b) added options to force gradient and chose the reorientation.   
#   c) opt- parameter is now read as whole line.   
#   d) added getXYZ to speed-up my current calculations.   
#   e) Re-added projector onto frequencies.   
#
#in version 0.1.5:  
#   a) Add function for reorientation of Cartesian Coordinates (Manipulate)      
#   b) removed the oldinvokelogging. I think, I will not need it (it was looking   
#       for whole as logging-level instead of just a letter)   
#   c) Manipulate is changed and still not working.   
#   d) Fixed some issues when working with gradients.   
#   e) changed Manipulate to MOI_reorient and changed the function.   
#
#in version 0.1.0:  
#   a) added some documentation/doc-strings  
#   b) added the functions printMat and printVec, taken from Duschinksy().  
#   c) fixed some errors  
#   d) made vector-output more beautiful  
#   e) try an oher invokeLogging; hopefully more stable.  
#   f) added warnings  

class  Spect:
   """ The class Spect is the parent-class for all spectrum-models in use here.
      Its init-function calculates/initialises most of the required data (the rest will
      be added in the derived-classes init-functions.)
      Some of the functions are redefined in the respective sub-classes.
         
      In addition, this function provides interfaces to the reading-routines (class Read)
      and the routines for generalising the (FC-based) spectra from One-Particle-Approx.
      to arbitrarily full spectrum.
      
      **USER-RELEVANT METHODS:**
      (There are many functions in this class that are intended to be called only by init.
      These classes will be changed to private and/or changed to be subclasses of __init__).
      The other functions are:  
      
      __init__(f)  ->(of course!)
         Its parameter is the input-file to Visper. From this, all other data are infered.
         Its returned variables are none.  

      calcspect()
         This function performs the actual calculation of the spectrum in the respect.
         model. This function is not present here at all!!! It is defined only for the
         derived classes (is there some virtual function in python!?)  
      
      outspect()
         This function performs the broadening of the line-spectrum calculated in 
         calcspect().  
      
      finish()
         This function maybe will be changed to be a destructor, if there exists such 
         in python. It cleans up everything uncleaned.
   
      **USER-RELEVANT MEMBERS**
      In principle, no relevant data should be required from here. But of course,
      many quantities are members of this class. These are, besides conversion factors,
      the energy, geometry,frequencies,...
   """
   # Below are the conversion factors and fundamental constant
   AMU2au=1822.88839                                          
   Angs2Bohr=1/0.52917721092                                  
   Hartree2GHz=6.579684e6                                     
   Hartree2cm_1=219474.63 
   
   # BEGIN OF DATA-DEF.
   Energy=[]
   dim=0
   CartCoord=[]
   read=[]  # this object is valid only after __init__. Maybe it doesn't belong here...
   spect=[[0],[0],[0]]
   width=5  # determines, how many rows of a matrix are writen aside each other.
   type="Spect"
   
   # The following additional elements are members of all instances:
   #
   #  logging 
   #  f  (frequencies in a 2-D array) #  Lmassw
   #  T  (temperature of the system)
   #  broadopt
   #  states1, states2
   #  read 
   #  mass (sqrt(masses))
   #  logfile
   #  F  mass-weighted hessian
 
   # END OF DATA-DEF.

   # BEGIN OF FUNCTION DEFINITIONS

   def __init__(self, f):
      """ This function initialises the Spect-object and creates its first.
         The INPUT, 'f', is the input-file to smallscript. In this function,
         only the information for output is read and the logging object,
         which has the output-file and the level of output-information is defined.
         All variables/functions that are common to all spectral tyes are initialised here.
      """
      def invokeLogging(logfile, mode="important"):
         """ initialises the logging-functionality
            **PARAMETERS**
            logfile   name of file to be used as log-file. It is expected to be an array of 
                   length 0 or one.
            mode:     5 different values are possible (see below); the lowest means: print 
                     much, higher values mean 
                     less printing

            **RETURNS**
            log:      opened file to write in
         """
         if logfile==[]:
            log=open("calculation.log", "a")
            self.logfile="calculation.log"
         else:
            s=logfile[-1].strip()
            log=open(s, "a")
            self.logfile=s
         
         #remove all white spaces and take only first character.
         mode= mode.strip()[0]
         #search, which option was set.
         if mode in ['a', "A", "0", 0]:
            logging=0
            log.write('use log-level all\n')
         elif mode in ['d', 'D', "1"]:
            logging=1
            log.write('use log-level detailed\n')
         elif mode in ['m', 'M','2']:
            logging=2
            log.write('use log-level medium\n')
         elif mode in ['i', 'I','3']:
            logging=3
         elif mode in ['s', 'S', '4']:
            logging=4
         else:
            logging=3
            log.write("logging-mode "+mode+" not recognized. Using 'important' instead\n")
         return logging, log
      
      #START LOG-FILE DEFINITION     
      #invoke logging (write results into file specified by 'out: ' or into 'calculation.log')
      logfile=re.findall(r"(?<=out: )[\w.,\(\) \=;:\-_]+", f, re.I)
      try:
         loglevel=re.findall(r"(?<=print:)[\w \d]+",f, re.I)[-1]
         self.logging = invokeLogging(logfile, loglevel )
      except IndexError:
         #else, try an other syntax
         try:
            loglevel=re.findall(r"(?<=print=)[\w \d]+",f, re.I)[-1]
            self.logging = invokeLogging(logfile, loglevel )
         except IndexError:
            self.logging = invokeLogging(logfile)
      
      # now,  write the header to the output-file.
      self.logging[1].write("\n==================================================================\n"
               "=====================  output of Visper  ========================\n"
               "==================================================================\n\n")
      self.logging[1].write("   INPUT-FILE:\n")
      self.logging[1].write(f)
      self.logging[1].write(" \n   END OF INPUT-FILE \n\n")
      #self.logging[1].write("calculations to be done: %s\n"%(todo))
      #END LOG-FILE DEFINITION     
      
      self.Energy=np.zeros(2)
      #START READING DATA FROM FILE
      # get files with electronic states:
      final=re.findall(r"(?<=final: )[\w.\-]+",f, re.I)
      initial=re.findall(r"(?<=initial: )[\w.\-]+",f, re.I)
      assert len(initial)==1, 'there must be one initial state'
      assert len(final)==1, 'there must be one final state'
      initial=initial[0]
      final=final[0]
      #check, if they are valid files and through an error if not.
      assert os.path.isfile(initial) and os.access(initial, os.R_OK),\
               initial+' is not a valid file name or not readable.'
      assert os.path.isfile(final) and os.access(final, os.R_OK),\
               final+' is not a valid file name or not readable.'
      
      #read options from input-file:
      opt=re.findall(r"(?<=opt:).*(?=\n)", f, re.M)
      if opt!=[]:
         self.opt=opt[-1]
      else:
         self.opt=" " 
         #make some empty string that the search-commands don't fail in ReadData.

      #find the information for the broadening
      broadopt=re.findall(r"(?<=broaden:)[\d\s\w.,\(\) \=;:\-_]+", f, re.M)
      if broadopt!=[]:
         self.broadopt=broadopt[-1]
      else:
         self.broadopt=" " 

      #find information on the systems temperature
      self.T=re.findall(r"(?<=T=)[\d .]+", self.opt, re.M)
      if self.T==[]:
         self.T=300
      else:
         self.T=float(self.T[-1])
      self.T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K
      
      if (re.search(r"(?<=width=)[\d .]+", self.opt, re.M)) is not None:
         self.width=int(re.findall(r"(?<=width=)[\d .]+", self.opt, re.M)[-1])

      # get number of maximal excited/ground from opt:
      states=re.findall(r"(?<=states=)[\d ,]*", self.opt, re.I)
      if len(states)==0:
         # i.e. no states given
	 self.states1=5
         self.states2=0
      else:
	 try:
            #i.e. 1 number given (else, the 'int' will give an error)
	    self.states1=int(states[0])
	    self.states2=0
	    self.logging[1].write("number of states: %d and %d \n"%(self.states1, self.states2))
	 except ValueError:
            try:
               # two numbers given, try to extract them
               self.states1=int(states[0].split(",")[0])
               self.states2=int(states[0].split(",")[1])
               self.logging[1].write("number of states: %d and %d \n"%(self.states1, self.states2))
            except ValueError:
               #unknown format. Use default and give a respective message.
               self.states1=5
               self.states2=0
               self.logging[1].write("!!number of vibrational states {0} is not an integer.",
                                    " Use default instead.\n".format(self.states1, self.states2))
      
      self.reader=Read.Read(initial, final)
      #The ReadData class finds out which format the input-files have
      # and reads the most important data from them.
      self.ReadData()

      #changes the orientation of the molecules to coincide with
      # the main axes of inertia. Input-orientations may
      # give too large or just unphysical HR-factors.
      manipulate=re.findall(r"(?<=manipulate:)[\w ,]*", self.opt, re.M)
      if manipulate!=[]:
         manipulate=manipulate[0].strip()
         if manipulate in ["moi", "MOI"]:
            self.MOI_reorient()
         if manipulate in ["rmsd" ,"RMSD"]:
            self.RMSD_reorient()
   
      #Calculate Frequencies and normal modes
      self.GetL()
      self.Duschinsky()
      
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
               #read file in format of spect
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
   
         spectfile=re.findall(r"(?<=spectfile=)[\w._\-]+", opt, re.I)
         if spectfile==[]:
            spectfile=re.findall(r"(?<=spectfile= )[\w._\-]+", opt, re.I)
            if spectfile==[]:
               spectfile=None
            else:
               spectfile=spectfile[-1]
         else:
            spectfile=spectfile[-1]
         return omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape, stick

      minint=0
      self.logging[1].write("\n STARTING TO CALCULATE BROADENED SPECTRUM.\n")
      omega, spectfile, gamma, gridpt, minfreq, maxfreq, shape, stick=handel_input(self.broadopt)
      #read file in format of spect
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
      if self.logging[0]<3:
         self.logging[1].write('neglect '+repr(minint)+' transitions, use only '+
                                repr(len(self.spect[1])-minint)+" instead.\n")
   
         if self.logging[0]<2:
            self.logging[1].write('minimal and maximal intensities:\n'+
              repr(self.spect[1][minint])+' '+repr(self.spect[1][-1])+"\n")
      
      #important for later loops: avoiding '.'s speeds python-codes up!!
      logwrite=self.logging[1].write  
     
      #make nPA from OPA if requested.
      n=re.findall(r"(?<=to nPA:)[ \d]*", self.broadopt, re.I)
      if n!=[]:
         MakeFull=MultiPart.OPAtoNPA(float(n[-1].strip()))
         self.logging[1].write("\n REACHING OPA TO NPA-PART. \n")
         self.logging[1].write(" ----------------------------------------"+
                                                 "-------- \n")
         MakeFull.GetSpect(self.spect, minint)
         TPAintens, TPAfreq=MakeFull.Calc()
      else: 
         TPAfreq=self.spect[0][minint:]
         TPAintens=self.spect[1][minint:]
      
      if stick:
         self.logging[1].write(" Intensity  \t frequency \n")
         for i in range(len(TPAfreq)):
            self.logging[1].write(" %3.6g  \t %3.6f\n"%(TPAintens[i],TPAfreq[i]))
            
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
      if self.logging[0]<3:
         self.logging[1].write('maximal and minimal frequencies:\n {0} {1}\n'.format(maxfreq, minfreq))
      #if no other grid is defined: use linspace in range
      if omega==None:
         omega=np.linspace(minfreq,maxfreq,gridpt)
         if self.logging[0]<2:
            self.logging[1].write("omega is equally spaced\n")
   
      sigma=gamma*2/2.355 #if gaussian used: same FWHM
      
      if gamma*1.1<=(maxfreq-minfreq)/gridpt:
         self.logging[1].write("\n WARNING: THE GRID SPACING IS LARGE COMPARED TO THE WIDTH OF THE PEAKS.\n"
              "THIS CAN ALTER THE RATIO BETWEEN PEAKS IN THE BROADENED SPECTRUM!")
   
      index=np.argsort(TPAfreq,kind='heapsort') #sort by freq
      freq=TPAfreq[index]
      intens=TPAintens[index]
   
      mini=0
      if spectfile==None:
         out=self.logging[1]
      else:
         out = open(spectfile, "w")
   
      if spectfile==None: #that means spectrum is printed into log-file
         logwrite("broadened spectrum:\n frequency      intensity\n")
      outwrite=out.write
      #this shrinks the size of the spectral lines; hopefully accelerates the script.
      #intens, freq=concise(intens,freq, sigma)
      lenfreq=len(freq)
      maxi=lenfreq-1 #just in case Gamma is too big or frequency-range too low
      for i in range(0,lenfreq-1):
         if freq[i]>=10*gamma+freq[0]:
            maxi=i
            break
      if shape=='g':
         sigmasigma=2.*sigma*sigma # these two lines are to avoid multiple calculations of the same
         npexp=np.exp
         intens/=sigma # scale all transitions according to their width.
         for i in xrange(len(omega)): 
            omegai=omega[i]
            for j in range(maxi,lenfreq):
               if freq[j]>=10*gamma+omegai:
                  maxi=j
                  break
            for j in range(mini,maxi):
               if freq[j]>=omegai-10*gamma:
                  # else it becomes -1 and hence the spectrum is wrong
                  mini=max(j-1,0) 
                  break
            spect=0
            for k in range(mini,maxi+1):
               spect+=intens[k]*npexp(-(omegai-freq[k])*(omegai-freq[k])/(sigmasigma))
            outwrite(u" %f  %e\n" %(omegai, spect))
      else:  #shape=='l':
         gammagamma=gamma*gamma
         for i in xrange(len(omega)): 
            omegai=omega[i]
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
            for k in range(mini,maxi+1):
               spect+=intens[k]*gamma/((omegai-freq[k])*(omegai-freq[k])+gammagamma)
            outwrite(u" %f   %e\n" %(omegai, spect))
      if spectfile!=None:
         #only close file if it was opened here
         out.close()

   def quantity(self, dim):
      """ Here some frequencies are defined; it is just for clearer code.
         This function is called by CalculationHR.
   
         **PARAMETERS**
         dim      dimension of matrices/vectors
      """
      F=np.zeros((2, dim, dim)) 
      self.CartCoord=np.zeros((2, 3, dim//3))
      P=np.zeros((2, dim,dim))
      return F, P
   
   def finish(self):
      """This simple function has the only task to close the log-file after 
         calculation as it is nice. 
         More over, it gives the user feedback about successfull finish;
         also nice if some warnings appeared before.
      """
      self.logging[1].close()
      #count warnings in self.logging[1]:
      logf=open(self.logfile,"r")
      warnings=0
      for line in logf:
         if 'WARNING:' in line:
            warnings+=1
      logf.close()
      foo=open(self.logfile, 'a')
      if warnings==0:
         foo.write("\n==================================================================\n"
                  "=========== VISPER FINISHED OPERATION SUCCESSFULLY.  =============\n"
                  "==================================================================\n\n")
      else:
         foo.write("\n==================================================================\n"
                  "============== VISPER FINISHED WITH "+repr(warnings)+" WARNINGS.  ===============\n"
                  "==================================================================\n\n")
      foo.close()
   # END OF FUNCTION DEFINITIONS
   
   def printMat(self,mat):
      """Function to print matrices in a nice way to the log-file.
         To keep it general, the matrix needs to be given as argument.
      """
      k=range(0,len(mat[0]))
      s=0
      t=min(s+self.width,len(mat[0]))
      while s<len(mat[0]):
         self.logging[1].write("\n")
         for temp in range(s,t):
            self.logging[1].write("          %03d "%(k[temp]+1))
         self.logging[1].write("\n")
         for j in range(len(mat)):
            self.logging[1].write(" %03d"%(j+1))
            for temp in range(s,t):
               self.logging[1].write("  %+.5e"%(mat[j][k[temp]]))
            self.logging[1].write("\n")
         s=t
         t=min(s+self.width,len(mat[0]))

   def printVec(self,vec):
      """This funcion is no that tricky but better than rewriting it everywhere it is
      indeed.
      """
      num = self.width//2
      self.logging[1].write("\n")
      if len(vec)>num:
         for j in range(len(vec)//num):
            for k in range(num):
               self.logging[1].write("    %03d  %e \t"%(j+k*len(vec)//num+1, vec[j+k*len(vec)//num]))
            self.logging[1].write("\n")
      else:
         for k in range(len(vec)):
            self.logging[1].write("    %03d  %e \t"%(k+1, vec[k]))
      self.logging[1].write("\n")

   def ReadData(self):
      """ This function gathers most essential parts for calculation of
         HR-factors from g09-files. That is: read neccecary information from the
         g09-files and calculate HR-factors as well as the  Duschinsky-rotation
         matrix and the shift between minima (needed later if the option Duschinsky
         is specified)
      """ 
      def IsZero(mat):
         """checks, if a matrix is zero
         """
         return np.all(mat==0)

      #read coordinates, force constant, binding energies from log-files and 
      # from the file, using the type of file that is known now...
      self.Energy[0], self.Energy[1]=self.reader.Energy()
      
      self.mass=self.reader.mass()
      if np.any(self.mass[0]==0):
         if np.any(self.mass[1]==0):
            assert 1==2, "Your input-files don't contain any mass."
         self.mass=self.mass[1]
      else:
         self.mass=self.mass[0]
      
      self.dim=len(self.mass)*3
      if self.logging[0]<2:
         self.logging[1].write("Dimensions: %d\n" %self.dim)
         self.logging[1].write(" Masses: \n")
         self.printVec(self.mass*self.mass) # print not the square-root!
   
      self.F=self.reader.Force()
      if IsZero(self.F[0]):
         self.F[0]=self.F[1]
         self.logging[1].write("WARNING: Only one force constant matrix given.\n")
         assert (self.type=='FC'), "You must specify both force constant matrices in this model."
      elif IsZero(self.F[1]):
         self.F[1]=self.F[0]
         self.logging[1].write("WARNING: Only one force constant matrix given.\n")
         assert (self.type=='FC'), "You must specify both force constant matrices in this model."
      assert not IsZero(self.F[0]), "ERROR: No force constant matrix given."
      
      self.CartCoord=self.reader.Coordinates()
      #this ugly statement is true if the molecule is just shifted (not rotated or deformed)
      # between states. I need to allow for shifts as it seems from looking at 
      # results obtained from NWChem. 
      # Moreover, I allow for very small changes that could occur due to shifting or
      # roundoff e.g. when converting units.
      if all(abs(self.CartCoord[0][i][j]-self.CartCoord[1][i][j] - 
                  self.CartCoord[0][i][0]+self.CartCoord[1][i][0]) <0.00001
                      for i in range(3) for j in range(self.dim//3) ):
         self.Grad=self.reader.Gradient() 
      elif re.search(r"gradient", self.opt, re.M) is not None:
         self.Grad=self.reader.Gradient() 
      else: 
         #define this for easier syntax in function MOI_reorient
         # and Duschinsky.
         self.Grad=[0,0]
      if self.logging[0]<3:
         self.logging[1].write('difference of minimum energy between states:'
                        ' Delta E= {0}\n'.format((self.Energy[0]-self.Energy[1])*self.Hartree2cm_1))
         if self.logging[0]<2:
            self.logging[1].write('Cartesian coordinates of initial state: \n')
            self.printMat(self.CartCoord[0].T/self.Angs2Bohr)
            self.logging[1].write('Cartesian coordinates of final state: \n')
            self.printMat(self.CartCoord[1].T/self.Angs2Bohr)
            if self.logging==0:
               self.logging[1].write('initial state: \n')
               self.printMat(F[0])
               self.logging[1].write('final state: \n')
               self.printMat(F[1])
   
   def makeXYZ(self):
      """This function creates an .xyz-file for opening e.g. with Chemcraft.
         It is intended to be called whenever it is likely that the input-states do not
         fit together or harmonic approximation will not hold.
         Hence, this is for backup only but in these cases might be helpful.
      """
      output=open("molecule.xyz", "w")

      #first, print the initial state:
      output.write("%d\n    inital state\n"%(self.dim//3))
      for i in range(len(self.CartCoord[0][0])):
         output.write("%d    %f   %f   %f\n"%(round(self.mass[i]*self.mass[i]/self.AMU2au/2),
                                             self.CartCoord[0][0][i]/self.Angs2Bohr,
                                             self.CartCoord[0][1][i]/self.Angs2Bohr,
                                             self.CartCoord[0][2][i]/self.Angs2Bohr) )
      #second, print the final state:
      output.write("\n%d\n    final state\n"%(self.dim//3))
      for i in range(len(self.CartCoord[0][0])):
         output.write("%d    %f   %f   %f\n"%(round(self.mass[i]*self.mass[i]/self.AMU2au/2),
                                             self.CartCoord[1][0][i]/self.Angs2Bohr,
                                             self.CartCoord[1][1][i]/self.Angs2Bohr,
                                             self.CartCoord[1][2][i]/self.Angs2Bohr) )

      #finally, print block with both states after another:
      output.write("\n%d\n    Both states\n"%(self.dim//3)*2)
      for i in range(len(self.CartCoord[0][0])):
         output.write("%d    %f   %f   %f\n"%(round(self.mass[i]*self.mass[i]/self.AMU2au/2),
                                             self.CartCoord[0][0][i]/self.Angs2Bohr,
                                             self.CartCoord[0][1][i]/self.Angs2Bohr,
                                             self.CartCoord[0][2][i]/self.Angs2Bohr) )
      for i in range(len(self.CartCoord[0][0])):
         output.write("%d    %f   %f   %f\n"%(round(self.mass[i]*self.mass[i]/self.AMU2au/2),
                                             self.CartCoord[1][0][i]/self.Angs2Bohr,
                                             self.CartCoord[1][1][i]/self.Angs2Bohr,
                                             self.CartCoord[1][2][i]/self.Angs2Bohr) )
      output.close()

   def GetL(self):
      """Function that calculates the frequencies and normal modes from force constant matrix.It
         mainly calls a hermitian eigenvalue-solver and computes the frequencies as 
         srqrt(eigenvalue) and the normal modes by cutting the translation/rotation off.
         A projection method once has been in use but turned out to be pointless.
   
         This function is a member of Spect but syntactically not bound to it at the moment.
         Hence, it has access to member-variables.
      
         **PARAMETERS** 
         self -> it is member of class
         F    -> two matrices, F[0] and F[1] ; the respective nuclear energy Hessians.
   
         **RETURN**
         f  -> frequencies of initial (f[0]) and final (f[1]) state. The dimesion is (2, dim-6)
         L  -> unitary matrices (L[0], L[1]) that diagonalise M*F (massweighted Hessian). Its columns are normal modes 
               in Cartesian Coordinats The dimension is (2, dim, dim-6)
         Lmass -> L*M where M_ij=1/sqrt(m_i m_j). This is mass-weighted L and will be used later for most systems.
        
      """
      # Defining arrays
      lenF=len(self.F[0])
      L=np.zeros(( 2, lenF, lenF-6 )) 
      Lmass=np.zeros(( 2, lenF, lenF-6 ))
      f=np.zeros(( 2, lenF-6 ))

      def gs(A):
         """This function does row-wise Gram-Schmidt orthonormalization of matrices. 
            code for Gram-Schmidt adapted from iizukak, see https://gist.github.com/iizukak/1287876
         """
         X=A.T # I want to orthogonalize row-wise
         Y = []
         npdot=np.dot
         for i in range(len(X)):
            temp_vec = X[i]
            for inY in Y :
               #proj_vec = proj(inY, X[i])
               proj_vec = map(lambda x : x *(npdot(X[i],inY) / npdot(inY, inY)) , inY)
               temp_vec = map(lambda x, y : x - y, temp_vec, proj_vec)
            Y.append( temp_vec/np.linalg.norm(temp_vec)) # normalise vectors
         return np.matrix(Y).T # undo transposition in the beginning

      def GetProjector(self, i):
         """ 
         """
         # Getting tensor of inertia, transforming to principlas axes
         moi=np.zeros((3,3))# this is Moment Of Inertia
         # print Coord[1]
         for j in [0,1,2]:
            for k in [0,1,2]:
               if k == j:
                  moi[j][j]=np.sum(self.mass*self.mass*(self.CartCoord[i][0]*self.CartCoord[i][0]+\
                           self.CartCoord[i][1]*self.CartCoord[i][1]+self.CartCoord[i][2]*self.CartCoord[i][2]-\
                           self.CartCoord[i][j]*self.CartCoord[i][k]))
               else:
                  moi[j][k]=np.sum(self.mass*self.mass*(self.CartCoord[i][j]*self.CartCoord[i][k]))
         diagI,Moi_trafo=np.linalg.eig(moi) # this can be shortened of course!
         index=np.argsort(diagI,kind='heapsort')
         #X=np.matrix(X[index]) #sorting by eigenvalues
         Moi_trafo=np.matrix(Moi_trafo) #notsorting by eigenvalues
      
         #now, construct the projector onto frame of rotation and translation using Sayvetz conditions
         D=np.zeros((self.dim,6))
         for k in [0,1,2]:# first three rows in D: The translational vectors
            for j in range(self.dim//3):
               #translations in mass-weighted coordinates
               D[3*j+k][k]=self.mass[j]
         for k in range(self.dim):# next three rows in D: The rotational vectors
            #rotations in mass-weighted coordinates
            D[k][3:6]=(np.cross(np.dot(Moi_trafo,self.CartCoord[i])[:].T[k//3],Moi_trafo[:].T[k%3]))*self.mass[k//3]
         D_orthog=gs(np.array(D)) #orhogonalize it
         ones=np.identity(self.dim)
         one_P=ones-np.dot(D_orthog,D_orthog.T)
         prob_vec=(D_orthog.T[1]+D_orthog.T[4]+D_orthog.T[0]+D_orthog.T[5]).T #what is this actually??
         assert not np.any(np.abs(prob_vec-np.dot(np.dot(D_orthog,D_orthog.T),prob_vec))>0.00001), \
                  'Translations and rotations are affected by projection operator.'+\
                  repr(np.abs(prob_vec-np.dot(np.dot(D_orthog,D_orthog.T),prob_vec)))
         assert not  np.any(np.abs(np.dot(one_P,prob_vec))>0.00001), \
                  "Projecting out translations and rotations from probe vector"
         return one_P
      
      def SortL(J,L,f):
         """This functions resorts the normal modes (L, f) such that the Duschinsky-Rotation
         matrix J becomes most close to unity (as possible just by sorting).
         In many cases, here chosing max(J[i]) does not help since there will be rows/columns occur
         more often. 
         Since I don't know any closed theory for calculating this, it is done by cosidering all possible cases.
         """
         
         #initialize the matrix that resorts the states:
         resort=np.zeros(np.shape(J))
         #FIRST, DO SOME GUESS HOW THE MATRIX COULD LOOK LIKE
         #chose largest elements in lines
         for i in range(len(J)):
            j=np.argmax(J[i])
            k=np.argmin(J[i])
            if J[i][j]>-J[i][k]:
               resort[i][j]=1
            else:
               resort[i][k]=1
         #now, go through rows and check if they are ok:
         #print "resortJ\n",resort
         resort=resort.T
         Nos=[]
         freePlaces=[]

         # NOW LOOK FOR ERRORS IN THE GUESS
         for i in range(len(J)):
            if sum(resort[i])==1:
               #this is the normal case: the order of
               # states did not change.
               continue
            elif sum(resort[i])==0:
               Nos.append(i)
            else:
               index=np.where(resort[i]==1)
               x=np.argmax(np.abs(J[index,i]))
               index=np.delete(index,x)
               resort[i][index]=0 #only x remains
               freePlaces=np.append(freePlaces,index)
         # By construction, this should always be true!
         assert len(Nos)==len(freePlaces), "dodododo!"
         freePlaces=np.array(freePlaces,dtype=int)
         
         #FIXING THE ERRORS NOW:
         #fill remaining lines. Therefore, set that element to one
         # whose value is largest under all remaining ones.
         # This method is not fair since the first has most choise but should
         # be fair enough in most cases
         for i in range(len(Nos)):
               x=np.argmax(np.abs(J[freePlaces,Nos[i]]))
               resort[Nos[i],freePlaces[x]]=1 #only x remains
               freePlaces=np.delete(freePlaces,x) # does it work the way I want it to work?
         assert len(freePlaces)==0, "the matrix is not readily processed."
         #FIX DONE.
         
         #since resort is a permutation matrix, it is unitary. Using this:
         return f.dot(resort), L.dot(resort)
         #  END OF SortL
  
      # do the following for both states:
      for i in [0,1]:
         # solve the eigenvalue-equation for F:

         #REMOVE ROTATIONS AND TRANSLATIONS FROM HESSIAN.
         #project=False
         project=True
         #project out the rotations and vibrations and not just throw away the smallest 6 eigen values
         # and respective eigen modes.
         if project:
            D=GetProjector(self, i)
            #ftemp,Ltemp=np.linalg.eig(np.dot(np.dot(D.T,self.F[i]),D))
            ftemp,Ltemp,info=dsyev(np.dot(np.dot(D.T,self.F[i]),D))
            #ftemp,Ltemp,info=dsyev(self.F[i])
            
            #Why can I not construct D such that I don't need to throw away anything?

            index=np.argsort(np.real(ftemp),kind='heapsort') # ascending sorting f
            # and then cut off the 6 smallest values: rotations and vibrations.
            f[i]=np.real(ftemp[index]).T[:].T[6:].T
            L[i]=(Ltemp.T[index].T)[:].T[6:].T
         else: #disabled right now, maybe I will reenable it later on!?
            #ftemp,Ltemp=np.linalg.eigh(self.F[i])
            ftemp,Ltemp,info=dsyev(self.F[i])
            for j in range(0,self.dim):
               norm=np.sum(Ltemp.T[j].T.dot(Ltemp.T[j]))
               print norm
               if np.abs(norm)>1e-12:
                  Ltemp.T[j]/=np.sqrt(norm)

            #sort the results with increasing frequency (to make sure it really is)
            index=np.argsort(np.real(ftemp),kind='heapsort') # ascending sorting f
            # and then cut off the 6 smallest values: rotations and vibrations.
            f[i]=np.real(ftemp[index]).T[:].T[6:].T
            L[i]=(Ltemp.T[index].T)[:].T[6:].T
         #END REMOVING ROTATIONS AND TRANSLATIONS.
   
         #the frequencies are square root of the eigen values of F
         # Here, we need to take care for values <0 due to bad geometry
         for j in range(len(f[i])):
            f[i][j]=np.sign(f[i][j])*np.sqrt(np.abs(f[i][j]))

         #through a warning if imaginary frequencies occured and take their pos. values for that.
         if np.any(f[i]<0):
            self.logging[1].write('WARNING: imaginary frequencies occured. The absolute'
                    ' values are used in the following.\n{0}\n'.format(f[i]))
            f[i]=np.abs(f[i])

         #calculate the mass matrix, for Lmass
         M=np.eye(self.dim)
         for j in range(self.dim):
            M[j,j]/=self.mass[j//3]
         Lmass[i]=M.dot(L[i])
         #if log level=0 or 1:
         if self.logging[0]<2:
            if i==0:
               self.logging[1].write("Initial states:\n")
            else:
               self.logging[1].write("Final states:\n")
            self.logging[1].write("Frequencies (cm-1)\n")
            self.printVec(f[i]*self.Hartree2cm_1)
            #indeed, this matrix is not resorted yet and yes, the user is not informed
            #that there are manipulations done on it afterwards.
            self.logging[1].write("L-matrix \n")
            self.printMat(Lmass[i])
      
      # to account for the effect that the frequencies independently change between the states 
      #  and hence may change their order, the final states L and f are resorted via the
      #  largest overlap. When strong changes occur, there is no fair method. This one tries
      #  to be as fair as possible.
      f[1],L[1]=SortL(np.linalg.pinv(L[0]).dot(L[1]),L[1],f[1])

      #recalculate Lmass for final state.
      Lmass[1]=M.dot(L[1]) # recalculate Lmass!
      self.Lmassw=Lmass
      self.f=f
      self.L=L
   
   def Duschinsky(self):
      """This function calculates the shift between two electronic states 
         (whose geometry is known, see x) as well as the
         Duschinsky-rotation matrix.

         **PARAMETERS:**
         mass:    array of square-roots of nuclear masses (length: N)

         **RETURN:**
         J:    Duschinsky-rotation matrix
         K:    displacement-vector of energy-minima in normal coordinates
      """
      self.J=np.zeros((self.dim-6,self.dim-6))
      self.K=np.zeros(self.dim-6)
      M=np.zeros((self.dim,self.dim))
   
      for i in range(self.dim):
         M[i][i]=self.mass[i//3] #square root of inverse masses
      #J=np.dot(L[0].T, L[1])  # for Lsorted
      self.J=np.linalg.pinv(self.Lmassw[0]).dot(self.Lmassw[1]) # for Lmassw
   
      #print "J\n", J
      if any(self.Grad[i]>0 for i in range(len(self.Grad))):
         self.K=self.Grad.T.dot(self.Lmassw[0])
         #self.K=(np.linalg.pinv(self.Lmassw[0]).dot(self.Grad)).T  # w p Lmassw
         # scale consistently: Now it is really the shift in terms of normal modes
         self.K/=self.f[0]*self.f[0]#*np.sqrt(2)  #
         #self.K=self.K[0] # it is matrix and needs to be a vector...
         #FIXME: this seems to be inconsistent! may be matrix or vector...

      else:
         DeltaX=np.array(self.CartCoord[1]-self.CartCoord[0]).flatten('F')  # need initial - final here.
         if self.logging[0] <1:
            self.logging[1].write('changes of Cartesian coordinates:\n')
            self.printVec(DeltaX)
         self.K=np.linalg.pinv(self.Lmassw[0]).dot(DeltaX)  # w p Lmassw
      
      if self.logging[0]<2:
         # print the Duschinsky matrix in a nice format
         self.logging[1].write('Duschinsky rotation matrix:\n')
         self.printMat(self.J)
         self.logging[1].write('\nDuschinsky displacement vector:\n')
         self.printVec(self.K)
 
 # FIXME: new class/file for the manipulation of data: 
 #        there is the L-matrix (reordering) and coordinate-changes (2 methods)

   def MOI_reorient(self):
      """This function reorients the final state in space such that
         the moment of inertia frames coincide.
         The function will fail when the order of the moments changes
         or are close to degenerate and hence are rotated towards each other.
         FIXME: I want to add a threshold to make the programme decide itself,
         whether this method is applicable or not.
      """
      #FIRST STEP: move the center of mass (COM) to origin:
      COM=np.zeros(3)
      X=np.zeros( (2,3,3) )
      diagI=np.zeros( (2,3) )
      #do it for initial and final state:
      for i in [0,1]:
         #loop over coordinates:
         for j in [0,1,2]:
            COM[j]=np.sum(self.CartCoord[i][j]*self.mass)
            COM[j]/=np.sum(self.mass) 
         #now it is Cartesian center of mass
         if self.logging[0]<2:
            if i==0:
               self.logging[1].write("Center of mass (initial) coordinates (Bohr):\n")
            else:
               self.logging[1].write("Center of mass (final) coordinates (Bohr):\n")
            self.printVec(COM)
         for j in [0,1,2]:
            #displacement of molecule into center of mass:
            self.CartCoord[i][j]-=COM[j]
      
         #SECOND STEP: Calculate the inertia-system.
         #  (MOI: Moment Of Inertia)
         MOI=np.zeros((3,3))# this is Moment Of Inertia
         #loop over coordinates:
         for j in [0,1,2]:
            MOI[j][j]=np.sum(self.mass*self.mass*(self.CartCoord[i][0]*self.CartCoord[i][0]+\
                     self.CartCoord[i][1]*self.CartCoord[i][1]+self.CartCoord[i][2]*self.CartCoord[i][2]-\
                     self.CartCoord[i][j]*self.CartCoord[i][j]))
            for k in range(j):
               MOI[j][k]=np.sum(self.mass*self.mass*(self.CartCoord[i][j]*self.CartCoord[i][k]))
               MOI[k][j]=MOI[j][k]
         #calculate the eigen-system of MOI
         diagI[i],X[i]=np.linalg.eig(MOI) 
         # sort it such, that the rotation is minimal. Therefore, first check
         # that it is not too big; otherwise, this simple code might fail.
         index=np.argsort(diagI[i], kind='heapsort')
         X[i]=(X[i].T[index]).T
         diagI[i]=diagI[i][index]
      #end for i in [0,1].
      
      #output on the information gathered above
      if self.logging[0]==0:
         self.logging[1].write('\nRotational constants (GHz) in principle axes\n')
         self.logging[1].write('initial state: '+ repr(1/(2*diagI[0].T)*self.Hartree2GHz)+"\n")
         self.logging[1].write('final state: '+ repr(1/(2*diagI[1].T)*self.Hartree2GHz)+"\n")
         self.logging[1].write("Inertia system in Cartesion Coordinates of initial state:\n")
         self.printMat(X[0])
         self.logging[1].write("Inertia system in Cartesion Coordinates of final state:\n")
         self.printMat(X[1])
      
      #overlap of the moi-systems: gives the rotation of the respective frames
      # but is free in sign; correct this in the following:
      O=X[0].dot(X[1].T) 

      rmsdi=[]
      # now: test all combinations of signs. criterion is the least square of 
      # Cartesian coordinates.
      sign=np.eye(3)
      for i in range(13):
         #indeed, this gives all combinations 
         # after another...
         sign[int((i//3+i)%3)][int((i//3+i)%3)]*=-1
         if i in [4,7,8,10, 11]:
            #they give redundant results.
            continue
         U=sign.dot(O.T)
         print sign
         rmsdi.append(RMSD(self.CartCoord[0]-U.dot(self.CartCoord[1])))

         #print  RMSD(self.CartCoord[0]-self.CartCoord[1]), RMSD(self.CartCoord[0]-U.dot(self.CartCoord[1]))
      rmsd=RMSD(self.CartCoord[0]-self.CartCoord[1])
      rmsdi.append(rmsd)
      #reassign i 
      i=np.argmin(rmsdi)
      if i==len(rmsdi)-1:
         # no rotation into MOI-frame should be done because initial 
         # geometry is the best. This is the case especially, if there
         # are (almost) degenerate moments of intertia.
         return
      #recover the sing-combination with least squares:
      if i>=4:
         i+=1
         if i>=7:
            i+=1
            if i>=8:
               i+=1
               if i>=10:
                  i+=1
                  if i>=11:
                     i+=1
      sign=np.eye(3)
      for j in range(i):
         sign[int((j//3+j)%3)][int((j//3+j)%3)]*=-1
      U=sign.dot(O)
      #apply this combination to coordinates of final state
      self.CartCoord[1]=U.dot(self.CartCoord[1])
      # produce respective matrix to rotate Force-constant matrix as well
      Y=np.zeros( (self.dim,self.dim) )
      for j in range(self.dim//3):
         Y[3*j:3*j+3].T[3*j:3*j+3]=U
      # apply the respective rotation:
      self.F[1]=np.dot(Y.T.dot(self.F[1]),Y)
      if any(self.Grad[i]>0 for i in range(len(self.Grad))):
         self.Grad=Y.dot(self.Grad)
   
      # finally: print what is done.
      if self.logging[0]==0:
         self.logging[1].write("Rotation of final state:\n")
         self.printMat(U)
         self.logging[1].write("Coordinates after manipulation::\n")
         self.logging[1].write('Cartesian coordinates of initial state: \n')
         self.printMat(self.CartCoord[0].T/self.Angs2Bohr)
         self.logging[1].write('Cartesian coordinates of final state: \n')
         self.printMat(self.CartCoord[1].T/self.Angs2Bohr)
   
   def RMSD_reorient(self):
      """This function reorients the final state in space such that
         the moment of inertia of coordinates is minimized.
         I assume that there might be coordinate flips and/or switches
         but no strong rotations by ~40 degree
      """
      #FIRST STEP: move the center of mass (COM) to origin:
      COM=np.zeros(3)
      X=np.zeros( (2,3,3) )
      diagI=np.zeros( (2,3) )
      #do it for initial and final state:
      for i in [0,1]:
         #loop over coordinates:
         for j in [0,1,2]:
            COM[j]=np.sum(self.CartCoord[i][j]*self.mass)
            COM[j]/=np.sum(self.mass) 
         #now it is Cartesian center of mass
         if self.logging[0]<2:
            if i==0:
               self.logging[1].write("Center of mass (initial) coordinates (Bohr):\n")
            else:
               self.logging[1].write("Center of mass (final) coordinates (Bohr):\n")
            self.printVec(COM)
         for j in [0,1,2]:
            #displacement of molecule into center of mass:
            self.CartCoord[i][j]-=COM[j]
      
      #SECOND STEP: check for all combinations of axis-flips and coordinate-switches.
      # test all combinations of signs. and permutations:
      rotated=self.CartCoord[1]
      for j in range(6):
         #loop over all permutations:
         O=np.zeros( (3,3) )
         if j%2==0:
            #for cyclic permutations
            for i in range(3):
               #for cyclic ones, this works...
               O[(i+j//2)%3][(i-j//2)%3]=1
         else:
            #for anti-cyclic permutations
            for i in range(3):
               #for cyclic ones, this works...
               O[-((i+j//2)%3)-1][(i-j//2)%3]=1
         sign=np.eye(3)
         for i in range(13):
            #loop over all sing-combinations.
            #indeed, this gives all combinations 
            # after another...
            sign[int((i//3+i)%3)][int((i//3+i)%3)]*=-1
            if i in [4,7,8,10, 11]:
               #they give redundant results.
               continue
            U=sign.dot(O)
            #print np.shape(rotated) ,np.shape(self.CartCoord[0]), np.shape(U.dot(self.CartCoord[1]))
            if RMSD(self.CartCoord[0]-rotated) >RMSD(self.CartCoord[0]-U.dot(self.CartCoord[1])):
               #if the current combination is the best, save it.
               rotated=U.dot(self.CartCoord[1])
      #save the flipped/switched system in CartCoord:
      self.CartCoord[1]=rotated

      #THIRD STEP: follow some Monte Carlo scheme.
      # comparatively small angles. Do 200 steps.
      U=np.eye(3)
      for i in range(40):
         #chose some angle:
         #for every level: do 5 tests and the refine the search...
         for j in range(50):
            #a rotation-matrix with small rotations:
            theta=0.063*(random.random()-.5)/(i+1)
            eta=0.063*(random.random()-.5)/(i+1)
            nu=0.063*(random.random()-.5)/(i+1)
            R_x=np.matrix([ [1,0,0],
                  [0,np.cos(theta),-np.sin(theta)],
                  [0,np.sin(theta),np.cos(theta)] ], dtype=float)
            R_y=np.matrix([[np.cos(eta),0,-np.sin(eta)],
                  [0,1,0],
                  [np.sin(eta),0,np.cos(eta)]], dtype=float)
            R_z=np.matrix([[np.cos(nu),-np.sin(nu),0],
                  [np.sin(nu),np.cos(nu),0],
                  [0,0,1] ], dtype=float)
            R=R_x.dot(R_y).dot(R_z)
            test=R.A.dot(self.CartCoord[1])
            if RMSD(self.CartCoord[0]-self.CartCoord[1])> RMSD(self.CartCoord[0]-test):
               #if it gets better: apply the change.
               self.CartCoord[1]=test
               U=R.A.dot(U)
      
      #FOURTH STEP: apply the rotation.

      # produce respective matrix to rotate Force-constant matrix as well
      Y=np.zeros( (self.dim,self.dim) )
      for j in range(self.dim//3):
         Y[3*j:3*j+3].T[3*j:3*j+3]=U
      # apply the respective rotation:
      self.F[1]=np.dot(Y.T.dot(self.F[1]),Y)
      if any(self.Grad[i]>0 for i in range(len(self.Grad))):
         self.Grad=Y.dot(self.Grad)
   
      # finally: print what is done.
      if self.logging[0]==0:
         self.logging[1].write("Rotation of final state::\n")
         self.printMat(U)
         self.logging[1].write("Coordinates after manipulation::\n")
         self.logging[1].write('Cartesian coordinates of initial state: \n')
         self.printMat(self.CartCoord[0].T/self.Angs2Bohr)
         self.logging[1].write('Cartesian coordinates of final state: \n')
         self.printMat(self.CartCoord[1].T/self.Angs2Bohr)

def RMSD(Delta):
   """This function calculates the RMSD of a matrix (intended
      for self.CartCoords)
   """
   rmsd=0
   for i in range(len(Delta)):
      for j in range(len(Delta[0])):
         rmsd+=Delta[i][j]*Delta[i][j]
   return rmsd

#version=0.1.6   

# End of Spect.py
