#!/usr/bin/python2
# filename: Spect.py
import numpy as np
import re, mmap, os.path, math, sys
#include [[Read.py]]
import Read
#include [[file_handler.py]]
import file_handler as logger
#include [[normal_modes.py]]
import normal_modes as NormalModes
#include [[atoms_align.py]]
import atoms_align as AtAl

# CHANGELOG
# =========
#in version 0.2.0:  
#   a) removed function 'quantity()'; it was not in use any more.
#   b) removed the function IsZero() since it didn't do anything
#   c) changed format of self.CartCoord
#
#in version 0.1.7:  
#   a) fixed rmsd-reorient: forgot to reorder Forces, transform of 
#       gradient was wrong. Lead to ceveir errors.  
#   b) outsourced logging to file_handler.  
#   c) outsourced calculation of normal modes to to normal_modes.  
#   d) outsourced coordinate-transformations.
#   e) outsourced broadening of spectrum
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
   sameF=False #remembers, whether two different force constant matrices are in use.
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
      #START DEFINITION OF LOG-FUNCTION
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
      #END DEFINITION OF LOG-FUNCTION
      
      self.Energy=np.zeros(2)
      #START READING DATA FROM FILE
      # get files with electronic states:
      final=re.findall(r"(?<=final: )[\w.\-]+",f, re.I)
      initial=re.findall(r"(?<=initial: )[\w.\-]+",f, re.I)
      assert len(initial)==1,'there must be one initial state'
      assert len(final)==1,  'there must be one final state'
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
      
      #The ReadData class finds out which format the input-files have
      # and reads all important data from them. 
      self.reader=Read.Read(initial, final) #initialise class
      self.ReadData()                       #do the actual work

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
      # an object of the respective class and do the calculations there.
      # Except for the frequencies, the normal-mode based quantities are
      # always taken from self.nm.
      self.nm=NormalModes.NormalMode(self)
      self.nm.GetL()
      self.nm.Duschinsky()
      #give a copy of the frequencies to Spect because they are used more 
      # frequently.
      self.f=self.nm.f
      # this is how to use the function to print normal modes of final state
      #self.log.printNormalModes(self,1) # use 0 for initial state.
      
   def ReadData(self):
      """ This function gathers most essential parts for calculation of
         HR-factors from g09-files. That is: read neccecary information from the
         g09-files and calculate HR-factors as well as the  Duschinsky-rotation
         matrix and the shift between minima (needed later if the option Duschinsky
         is specified)
      """ 
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
      # print not the square-root! Print in amu for convenience.
      self.log.printVec(self.mass*self.mass/self.AMU2au)
   
      self.F=self.reader.Force()
      if np.all(self.F[0]==0):
         self.F[0]=self.F[1]
         self.log.write("WARNING: Only one force constant matrix given.\n")
         self.sameF=True
         assert (self.type=='FC'), "You must specify both force constant matrices in this model."
      elif np.all(self.F[1]==0):
         self.F[1]=self.F[0]
         self.log.write("WARNING: Only one force constant matrix given.\n")
         self.sameF=True
         assert (self.type=='FC'), "You must specify both force constant matrices in this model."
      assert not np.all(self.F[0]==0), "ERROR: No force constant matrix given."
      
      self.CartCoord=self.reader.Coordinates()
      #this ugly statement is true if the molecule is just shifted (not rotated or deformed)
      # between states. I need to allow for shifts as it seems from looking at 
      # results obtained from NWChem. 
      # Moreover, I allow for very small changes that could occur due to shifting or
      # roundoff e.g. when converting units.
      if all(abs(self.CartCoord[0][i*3+j]-self.CartCoord[1][i*3+j] - 
                  self.CartCoord[0][j]+self.CartCoord[1][j]) <0.00001
                      for j in range(3) for i in range(self.dim//3) ):
         self.Grad=self.reader.Gradient() 
      elif re.search(r"gradient", self.opt, re.M) is not None:
         self.Grad=self.reader.Gradient() 
         #self.log.write("gradient in input-format")
         #self.log.printVec(self.Grad)
      else: 
         #define this for easier syntax in function MOI_reorient
         # and Duschinsky.
         self.Grad=[0,0]
      if self.Energy[0]-self.Energy[1]<0:
         self.log.write('vertical relaxation energy:'
                        ' Delta E= {0}\n'.format((self.Energy[0]-self.Energy[1])*self.Hartree2cm_1), 3)
      else:
         self.log.write('vertical excitation energy:'
                        ' Delta E= {0}\n'.format((self.Energy[0]-self.Energy[1])*self.Hartree2cm_1), 3)
      if self.log.level<2:
            self.log.write('Cartesian coordinates of initial state: \n')
            self.log.printCoordinates(self.CartCoord[0]/self.Angs2Bohr)
            self.log.write('Cartesian coordinates of final state: \n')
            self.log.printCoordinates(self.CartCoord[1]/self.Angs2Bohr)
            if self.log.level==0:
               self.log.write('Hessian of initial state: \n')
               self.log.printMat(self.F[0])
               self.log.write('Hessian of final state: \n')
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
      for i in range(len(self.CartCoord[0]//3)):
         output.write("%d    %f   %f   %f\n"%(round(self.mass[i]*self.mass[i]/self.AMU2au/2),
                                             self.CartCoord[0][i*3+0]/self.Angs2Bohr,
                                             self.CartCoord[0][i*3+1]/self.Angs2Bohr,
                                             self.CartCoord[0][i*3+2]/self.Angs2Bohr) )
      #second, print the final state:
      output.write("\n%d\n    final state\n"%(self.dim//3))
      for i in range(len(self.CartCoord[0]//3)):
         output.write("%d    %f   %f   %f\n"%(round(self.mass[i]*self.mass[i]/self.AMU2au/2),
                                             self.CartCoord[1][i*3+0]/self.Angs2Bohr,
                                             self.CartCoord[1][i*3+1]/self.Angs2Bohr,
                                             self.CartCoord[1][i*3+2]/self.Angs2Bohr) )

      #finally, print block with both states after another:
      output.write("\n%d\n    Both states\n"%(self.dim//3)*2)
      for i in range(len(self.CartCoord[0]//3)):
         output.write("%d    %f   %f   %f\n"%(round(self.mass[i]*self.mass[i]/self.AMU2au/2),
                                             self.CartCoord[0][i*3+0]/self.Angs2Bohr,
                                             self.CartCoord[0][i*3+1]/self.Angs2Bohr,
                                             self.CartCoord[0][i*3+2]/self.Angs2Bohr) )
      for i in range(len(self.CartCoord[0]//3)):
         output.write("%d    %f   %f   %f\n"%(round(self.mass[i]*self.mass[i]/self.AMU2au/2),
                                             self.CartCoord[1][i*3+0]/self.Angs2Bohr,
                                             self.CartCoord[1][i*3+1]/self.Angs2Bohr,
                                             self.CartCoord[1][i*3+2]/self.Angs2Bohr) )
      output.close()
   # END OF FUNCTION DEFINITIONS

version='0.2.0'   
# End of Spect.py
