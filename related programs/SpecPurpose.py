#!/usr/bin/python2
# filename: SpecPurpose.py

from math import ceil
#include [[Spect.py]]:
import Spect

import numpy as np

# CHANGELOG
# =========
#in version 0.0.1:  
#   1) invented class VibDeform
#

class  VibDeform:
   def __init__(self, parent, state=0):
      """ Initialises the class:  
          Copies the required data from the Spect-class.
      """
      self.coords=parent.CartCoord[state]/parent.Angs2Bohr
      self.Modes=parent.nm.Lmassw[state] # it is not normalised!
      self.mass=parent.mass*parent.mass/parent.AMU2au
      #construct the ouptut-file:
      self.logwrite=parent.log.logfile.split(".")[0]

   def printDeformed(self, i, alpha):
      """Function to be called by the user: 
         ==PARAMETERS==
         i:      single integer, specifying the mode number to follow.
                 (later maybe will be generalised to support multi-dim.)
         alpha:  vector of floats, specifying the scaling of the mode
                 where the output should be given.
      """
      self.output=open(self.logwrite+".vibs", "a")
      for j,a in enumerate(alpha):
         self.__printvib(self.coords+a*self.Modes.T[i], a,i)
      self.output.close()

   def __printvib(self,geometry, deflection, mode_nr):
      """This function prints the normal modes (in Cartesian basis)
         to an extra file (*.nm) in the G09-log format. These files are understood
         by Chemcraft.
         ==PARAMETERS==
         geometry:   geometry to be printed
         deflection: the value of alpha, this is made with
         mode_nr:    the number of the mode that is changed.
      """

      # introducing the geometry:
      #GEOMETRY
      self.output.write("                Deflection: %f in mode nr. %i :\n"%(deflection, mode_nr))
      self.output.write(" ---------------------------------------------------------------------\n"+
        " Center     Atomic      Atomic             Coordinates (Angstroms)\n"+
        " Number     Number       Type             X           Y           Z\n"+
        " ---------------------------------------------------------------------\n")
      #print molecular geometry:
      for i in range(len(self.mass)):
         self.output.write("    %3d"%(i+1))
         self.output.write("        %3d"%ceil(self.mass[i]/2.)) 
         self.output.write("           0") 
         self.output.write("        %.6f    %.6f   %.6f\n"%(geometry[3*i+0],geometry[3*i+1],geometry[3*i+2]))
      self.output.write(" ---------------------------------------------------------------------\n\n")


version='0.0.1'   
# End of SpecPurpose.py
