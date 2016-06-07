#!/usr/bin/python2
# filename: file_handler.py
import math  # needed for function ceil in printNormalModes
from numpy import ceil

# CHANGELOG
# =========
#in version 0.1.0:  
#  1) Intialised class
#  2) Added function printNormalModes()
#  3) fixed error in printVec() to print not all values
#

class logging():
   """This class does everything related to output of the printed information.
      Its main job is to hold an open file handler (loghandler) and to write texts
      into a log-file.
   """
   logfile=''
   level='0'
   loghandler=""
   numWarn=0
   width=5  # determines, how many rows of a matrix are writen aside each other.
   # also, vectors are subdivided into width//2 blocks that are written aside each other.

   def __init__(self, level, logfile):
      """initialise the quantities important for writing.
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
      
      self.logfile=logfile
      self.level, self.loghandler = invokeLogging(logfile, level)
      self.write("\n==================================================================\n"
               "=====================  output of Visper  =========================\n"
               "==================================================================\n\n")

   def __del__(self):
      """This simple function has the only task to close the log-file after 
         calculation as it is nice. 
         More over, it gives the user feedback about successfull finish;
         also nice if some warnings appeared before.
      """
      self.loghandler.close()
      #count warnings in self.loghandler:
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
                  "=============== VISPER FINISHED WITH "+repr(warnings)+" WARNINGS.  ================\n"
                  "==================================================================\n\n")
      foo.close()
   
   def write(self, text, level=70):
      """ write a given text to the log-file if its priority (level) is larger
         than the threshold-priority for printing information
      """
      if level>self.level:
         self.loghandler.write(text)
   
   def printCoordinates(self,Coords):
      """This function prints the nuclear coordinates, given as a vector
         of the form
         (atom1_x, atom1_y, atom1_z, atom2_x,...,atomN_z)
         in a convenient format with the atoms as lines and coordinates
         as rows.
      """
      self.loghandler.write("\n")
      self.loghandler.write("          X               \t")
      self.loghandler.write("     Y               \t")
      self.loghandler.write("     Z\n")
      for j in range(len(Coords)//3):
         for k in [0,1,2]:
            self.loghandler.write(" %03d  %e \t"%(j+k*len(Coords)//3+1, Coords[j*3+k]))
         self.loghandler.write("\n")
      self.loghandler.write("\n")

   def printNormalModes(self,parent,i):
      """This function prints the normal modes (in Cartesian basis)
         to an extra file (*.nm) in the G09-log format. These files are understood
         by Chemcraft.
         ==PARAMETERS==
          parent - pointer to the Spect-class (or respectively inherited one)
          i      - state, for which the normal modes should be printed 
                   i=0: initial state
                   i=1: final state
      """

      #first, copy the quantities needed from the Spect-class
      # and normalise the mass-weighted transformation matrix locally.
      L=parent.nm.Lmassw[i].T
      coords=parent.CartCoord[i]/parent.Angs2Bohr
      mass=parent.mass*parent.mass/parent.AMU2au
      freq=parent.f[i]*parent.Hartree2cm_1
      #normalise the L-matrix:
      for i in range(len(L)):
         norm=np.sum(L[i]*L[i])
         L[i]=L[i]/np.sqrt(norm)
      L=L.T
      
      #construct the ouptut-file:
      nm_file=self.logfile.split(".")[0]
      if i==0:
         output=open(nm_file+"_init.nm", "w")
      else:
         output=open(nm_file+"_final.nm", "w")

      #print the header of the files for Chemcraft to recognise the G09-format
      #HEADER
      output.write(" Entering Gaussian System\n") 
      output.write(" ----------------------------------------------------------------------\n"+
           " #P BLYP/6-31G(d) Freq\n -----------------------------------------------------------------------\n")
      
      # introducing the geometry:
      #GEOMETRY
      output.write("                          Input orientation:\n")
      output.write(" ---------------------------------------------------------------------\n"+
        " Center     Atomic      Atomic             Coordinates (Angstroms)\n"+
        " Number     Number       Type             X           Y           Z\n"+
        " ---------------------------------------------------------------------\n")
      #print molecular geometry:
      for i in range(len(coords)):
         output.write("    %3d"%(i+1))
         output.write("        %3d"%math.ceil(mass[i]/2.)) 
         output.write("           0") 
         output.write("        %.6f    %.6f   %.6f\n"%(coords[3*i+0],coord[3*i+1],coord[3*i+2]))
      output.write(" ---------------------------------------------------------------------\n\n")
     
      #print the frequencies:
      #NORMAL MODES
      s=0
      t=min(s+3,len(L[0]//3))
      while s<len(L[0]//3):
         for k in range(s,t):
            output.write("                   %3d "%(k+1))
         output.write("\n")
         for k in range(s,t):
            output.write("                     A ")
         output.write("\n Frequencies --    ")
         for k in range(s,t):
            output.write("%3.4f               "%(freq[k]))
         output.write("\n")
         #some dummy-values that need to be there but are not used. Set them 0 therefore.
         output.write(" Red. masses --     0.0000                 0.0000                 0.0000\n")
         output.write(" Frc consts  --     0.0000                 0.0000                 0.0000\n")
         output.write(" IR Inten    --     0.0000                 0.0000                 0.0000\n")
         output.write("  Atom  AN")
         for k in range(s,t):
            output.write("      X      Y      Z  ")
         output.write("\n")
         for i in range(len(L)//3):
            output.write("   %3d"%(i+1))
            output.write(" %3d"%math.ceil(mass[i]/2.))
            for k in range(s,t):
               output.write("     %3.2f   %3.2f  %3.2f"%(L[3*i][k],L[3*i+1][k],L[3*i+2][k]))
            output.write("\n")
         s=t
         t=min(s+3,len(L[0]))
      output.close()

   def printVec(self,vec):
      """This funcion is no that tricky but better than rewriting it everywhere it is
         indeed.
      """
      num = self.width//2
      self.loghandler.write("\n")
      rows=int(ceil(float(len(vec))/num))
      if len(vec)>num:
         #for j in range(len(vec)//num):
         for j in range(rows):
            for k in range(num):
               if (j+k*rows)<len(vec):
                  self.loghandler.write("    %03d  %e \t"%(j+k*rows+1, vec[j+k*rows]))
               else:
                  #vector is out of range.
                  break
            self.loghandler.write("\n")
      else:
         for k in range(len(vec)):
            self.loghandler.write("    %03d  %e \t"%(k+1, vec[k]))
      self.loghandler.write("\n")

   def printMat(self,mat):
      """Function to print matrices in a nice way to the log-file.
         To keep it general, the matrix needs to be given as argument.
      """
      k=range(0,len(mat[0]))
      s=0
      t=min(s+self.width,len(mat[0]))
      while s<len(mat[0]):
         self.loghandler.write("\n")
         for temp in range(s,t):
            self.loghandler.write("          %03d "%(k[temp]+1))
         self.loghandler.write("\n")
         for j in range(len(mat)):
            self.loghandler.write(" %03d"%(j+1))
            for temp in range(s,t):
               self.loghandler.write("  %+.5e"%(mat[j][k[temp]]))
            self.loghandler.write("\n")
         s=t
         t=min(s+self.width,len(mat[0]))

version='0.1.0'
# End of file_handler.py
