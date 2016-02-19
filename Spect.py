#!/usr/bin/python2
# filename: Spect.py
import numpy as np
from scipy.linalg.lapack import dsyev
import re, mmap, os.path, math, sys
import Read
import MultiPart

Hartree2GHz=6.579684e6                                     
Hartree2cm_1=219474.63 

# CHANGELOG
# =========
#to version 0.1.5:  
#to version 0.1.0:  
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
   #  logging #  f  (frequencies in a 2-D array) #  Lmassw
   #  T  (temperature of the system)
   #  broadopt
   #  states1, states2
   #  read 
   #  mass (sqrt(masses))
   #  logfile
 
   # END OF DATA-DEF.

   # BEGIN OF FUNCTION DEFINITIONS

   def __init__(self, f):
      """ This function initialises the Spect-object and creates its first.
         The INPUT, 'f', is the input-file to smallscript. In this function,
         only the information for output is read and the logging object,
         which has the output-file and the level of output-information is defined.
         All variables/functions that are common to all spectral tyes are initialised here.
      """
      
      def oldinvokeLogging(logfile, mode="important"):
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
         else:
            s=logfile[-1].strip()
            log=open(s, "a")

         #remove all white spaces
         mode= mode.strip()
         #search, which option was set.
         if mode in ['all', "ALL", "All", "0", 0]:
            logging=0
            log.write('use log-level all\n')
         elif mode in ['detailed', 'DETAILED', 'Detailed', "1"]:
            logging=1
            log.write('use log-level detailed\n')
         elif mode in ['medium', 'MEDIUM','Medium', '2']:
            logging=2
            log.write('use log-level medium\n')
         elif mode in ['important', 'IMPORTANT','Important', '3']:
            logging=3
         elif mode in ['short', 'SHORT','Short', '4']:
            logging=4
         else:
            logging=3
            log.write("logging-mode "+mode+" not recognized. Using 'important' instead\n")
         print logging
         return logging, log
      
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
      opt=re.findall(r"(?<=opt:)[\d\s\w.,\(\) \=;:-]+", f, re.M)
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
               logging[1].write("number of states: %d and %d \n"%(states1, states2))
            except ValueError:
               #unknown format. Use default and give a respective message.
               self.states1=5
               self.states2=0
               logging[1].write("!!number of vibrational states {0} is not an integer.",
                                    " Use default instead.\n".format(self.states1, self.states2))
      
      self.reader=Read.Read(initial, final)
      #The ReadData class finds out which format the input-files have
      # and reads the most important data from them.
      self.ReadData()
      
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
      warnings=re.findall('WARNING:',self.logfile)
      warnings=len(warnings)
      foo=open(self.logfile, 'a')
      if warnings==0:
         foo.write("\n==================================================================\n"
                  "=========== VISPER FINISHED OPERATION SUCCESSFULLY.  =============\n"
                  "==================================================================\n\n")
      else:
         foo.write("\n==================================================================\n"
                  "========= VISPER FINISHED WITH "+repr(warnings)+" WARNINGS.  ===========\n"
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
      for j in range(len(vec)//num):
         for k in range(num):
            self.logging[1].write("    %03d  %e \t"%(j+k*len(vec)//num+1, vec[j+k*len(vec)//num]))
         self.logging[1].write("\n")
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
      assert IsZero(self.mass[0]-self.mass[1]), "It seems as if the masses don't coincide within files."
      self.mass=self.mass[0]
      
      self.dim=len(self.mass)*3
      if self.logging[0]<2:
         self.logging[1].write("Dimensions: %d\n" %self.dim)
         self.logging[1].write(" Masses: \n")
         self.printVec(self.mass)
   
      self.F=self.reader.Force()
      if IsZero(self.F[0]):
         self.F[0]=self.F[1]
         self.write("WARNING: Only one force constant matrix given.")
         assert (self.type=='FC'), "You must specify both force constant matrices in this model."
      elif IsZero(self.F[1]):
         self.F[1]=self.F[0]
         self.write("WARNING: Only one force constant matrix given.")
         assert (self.type=='FC'), "You must specify both force constant matrices in this model."
      assert not IsZero(self.F[0]), "ERROR: No force constant matrix given."
      
      self.CartCoord=self.reader.Coordinates()
      if IsZero(self.CartCoord[0]-self.CartCoord[1]):
         self.Grad=self.reader.Gradient() 
      
      if self.logging[0]<3:
         self.logging[1].write('difference of minimum energy between states:'
                        ' Delta E= {0}\n'.format((self.Energy[0]-self.Energy[1])*Hartree2cm_1))
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
   
      #Calculate Frequencies and normal modes
      self.GetL()
      self.Duschinsky()
   
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
      for i in range(2):
         # solve the eigenvalue-equation for F:

         # here one can choose between the methods: result is more or less 
         #  independent
         #ftemp,Ltemp=np.linalg.eig(self.F[i])
         #ftemp,Ltemp=np.linalg.eigh(self.F[i])
         ftemp,Ltemp,info=dsyev(self.F[i])
         
         #sort the results with increasing frequency (to make sure it really is)
         index=np.argsort(np.real(ftemp),kind='heapsort') # ascending sorting f
         # and then cut off the 6 smallest values: rotations and vibrations.
         f[i]=np.real(ftemp[index]).T[:].T[6:].T
         L[i]=(Ltemp.T[index].T)[:].T[6:].T
   
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
            self.printVec(f[i]*Hartree2cm_1)
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
      DeltaX=np.array(self.CartCoord[1]-self.CartCoord[0]).flatten('F')  # need initial - final here.
      if self.logging[0] <1:
         self.logging[1].write('changes of Cartesian coordinates:\n')
         self.printVec(DeltaX)
    #  K=(DeltaX.dot(M)).dot(L[0])
      #K=M.dot(L[0]).T.dot(DeltaX)  # with Lsorted
      self.K=np.linalg.pinv(self.Lmassw[0]).dot(DeltaX)  # w p Lmassw
      #K=L[0].T.dot(M).dot(M).dot(DeltaX)  # with Lmassw
   
      #K*=np.sqrt(np.pi)/2. #correction factor due to some magic reason  --> want to get rid of this!!!
   
      if self.logging[0]<2:
         # print the Duschinsky matrix in a nice format
         self.logging[1].write('Duschinsky rotation matrix:\n')
         self.printMat(self.J)
         self.logging[1].write('\nDuschinsky displacement vector:\n')
         self.printVec(self.K)

#version=0.1.5
# End of Spect.py
