#!/usr/bin/python2
# filename: Spect.py
import numpy as np
import re, mmap, os.path, math, sys
import Read
import MultiPart
import file_handler as logger
import normal_modes as NormalModes
import atoms_align as AtAl

# CHANGELOG
# =========
#in version 0.2.0:  
#
#
#in version 0.1.7:  
#   a) fixed rmsd-reorient: forgot to reorder Forces, transform of 
#       gradient was wrong. Lead to ceveir errors.  
#   b) outsourced logging to file_handler.  
#   c) outsourced calculation of normal modes to to normal_modes.  
#   c) outsourced coordinate-transformations.
#
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
   type="Spect"
   
   # The following additional elements are members of all instances:
   #
   #  f  (frequencies in a 2-D array) #  Lmassw
   #  T  (temperature of the system)
   #  broadopt
   #  states1, states2
   #  read 
   #  mass (sqrt(masses))
   #  log  handler of class logging
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
      #START LOG-FILE DEFINITION     
      #invoke logging (write results into file specified by 'out: ' or into 'calculation.log')
      logfile=re.findall(r"(?<=out: )[\w.,\(\) \=;:\-_]+", f, re.I)
      try:
         loglevel=re.findall(r"(?<=print:)[\w \d]+",f, re.I)[-1]
      except IndexError:
         #else, try an other syntax
         try:
            loglevel=re.findall(r"(?<=print=)[\w \d]+",f, re.I)[-1]
         except IndexError:
            loglevel='important'
      self.log=logger.logging(loglevel,logfile)
      
      # now,  write the header to the output-file.
      self.log.write("   INPUT-FILE:\n")
      self.log.write(f)
      self.log.write(" \n   END OF INPUT-FILE \n\n")
      #self.log.write("calculations to be done: %s\n"%(todo))
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
         self.log.width=int(re.findall(r"(?<=width=)[\d .]+", self.opt, re.M)[-1])

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
	    self.log.write("number of states: %d and %d \n"%(self.states1, self.states2))
	 except ValueError:
            try:
               # two numbers given, try to extract them
               self.states1=int(states[0].split(",")[0])
               self.states2=int(states[0].split(",")[1])
               self.log.write("number of states: %d and %d \n"%(self.states1, self.states2))
            except ValueError:
               #unknown format. Use default and give a respective message.
               self.states1=5
               self.states2=0
               self.log.write("!!number of vibrational states {0} is not an integer.",
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
         #initialise object of the respective class
         self.manipulate=AtAl.align_atoms(manipulate, self)
         #perform the calculation
         self.manipulate.perform()
      else:
         #initialise object of the respective class
         self.manipulate=AtAl.align_atoms(manipulate, self)
         # just shift the molecules both to their center of mass
         self.manipulate.shift()
         #copy the manipulated data back here.
         self.CartCoord=self.manipulate.CartCoord

   
      #Calculate Frequencies and normal modes. Therefore, intialise
      # an object of the respective class and do the respect. calculations
      self.nm=NormalModes.NormalMode(self.log, self.F, self.mass, self.CartCoord, self.Grad)
      self.nm.GetL()
      self.nm.Duschinsky()
      #give me a copy of the frequencies
      self.f=self.nm.f
      
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
      self.log.write("\n STARTING TO CALCULATE BROADENED SPECTRUM.\n")
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
      self.log.write('neglect '+repr(minint)+' transitions, use only '+
                                repr(len(self.spect[1])-minint)+" instead.\n", 3)
   
      self.log.write('minimal and maximal intensities:\n'+
              repr(self.spect[1][minint])+' '+repr(self.spect[1][-1])+"\n", 2)
      
      #important for later loops: avoiding '.'s speeds python-codes up!!
      logwrite=self.log.write  
     
      #make nPA from OPA if requested.
      n=re.findall(r"(?<=to nPA:)[ \d]*", self.broadopt, re.I)
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
      
      if stick:
         self.log.write(" Intensity  \t frequency \n")
         for i in range(len(TPAfreq)):
            self.log.write(" %3.6g  \t %3.6f\n"%(TPAintens[i],TPAfreq[i]))
            
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
      self.log.write('maximal and minimal frequencies:\n {0} {1}\n'.format(maxfreq, minfreq), 3)
      #if no other grid is defined: use linspace in range
      if omega==None:
         omega=np.linspace(minfreq,maxfreq,gridpt)
         self.log.write("omega is equally spaced\n",2)
   
      sigma=gamma*2/2.355 #if gaussian used: same FWHM
      
      if gamma*1.1<=(maxfreq-minfreq)/gridpt:
         self.log.write("\n WARNING: THE GRID SPACING IS LARGE COMPARED TO THE WIDTH OF THE PEAKS.\n"
              "THIS CAN ALTER THE RATIO BETWEEN PEAKS IN THE BROADENED SPECTRUM!")
   
      index=np.argsort(TPAfreq,kind='heapsort') #sort by freq
      freq=TPAfreq[index]
      intens=TPAintens[index]
   
      mini=0
      if spectfile==None:
         out=self.log
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
      self.log.write("Dimensions: %d\n" %self.dim,2)
      self.log.write(" Masses: \n",2)
      self.log.printVec(self.mass*self.mass) # print not the square-root!
   
      self.F=self.reader.Force()
      if IsZero(self.F[0]):
         self.F[0]=self.F[1]
         self.log.write("WARNING: Only one force constant matrix given.\n")
         assert (self.type=='FC'), "You must specify both force constant matrices in this model."
      elif IsZero(self.F[1]):
         self.F[1]=self.F[0]
         self.log.write("WARNING: Only one force constant matrix given.\n")
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
      self.log.write('difference of minimum energy between states:'
                        ' Delta E= {0}\n'.format((self.Energy[0]-self.Energy[1])*self.Hartree2cm_1), 3)
      if self.log.level<2:
            self.log.write('Cartesian coordinates of initial state: \n')
            self.log.printMat(self.CartCoord[0].T/self.Angs2Bohr)
            self.log.write('Cartesian coordinates of final state: \n')
            self.log.printMat(self.CartCoord[1].T/self.Angs2Bohr)
            if self.log.level==0:
               self.log.write('initial state: \n')
               self.log.printMat(self.F[0])
               self.log.write('final state: \n')
               self.log.printMat(self.F[1])
   
   def makeXYZ(self):
      """This function creates an .xyz-file for opening e.g. with Chemcraft.
         It is intended to be called whenever it is likely that the input-states do not
         fit together or harmonic approximation will not hold.
         Hence, this is for backup only but in these cases might be helpful.
      """
      molfile=self.log.logfile.split(".")[0]
      output=open(molfile+".xyz", "w")

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
 
   # END OF FUNCTION DEFINITIONS

#version=0.2.0   
# End of Spect.py
