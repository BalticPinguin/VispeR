#!/usr/bin/python2
# filename: HR_spects.py
import numpy as np
import re, mmap, os.path, math, sys
#include [[file_handler.py]]
import file_handler as logger

#include [[FC_spects.py]]
import FC_spects

# CHANGELOG
# =========
#in version 0.0.1:  
#  1) initialisation of class

class HR_spect(FC_spects.FC_spect):
   """ This class is intended to calculate FC-spectra using HR-factors and frequencies
       read from a file.
   """
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
      hr_file=re.findall(r"(?<=HR-file: )[\w.\-]+",f, re.I)
      assert len(hr_file)==1,'there must be one initial state'
      self.hr_file=hr_file[0]
      #check, if they are valid files and through an error if not.
      assert os.path.isfile(self.hr_file) and os.access(self.hr_file, os.R_OK),\
               self.hr_file+' is not a valid file name or not readable.'
      
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
      
         #initialise object of the respective class
         self.manipulate=AtAl.align_atoms(manipulate, self)
         # just shift the molecules both to their center of mass
         self.manipulate.shift()
         #copy the manipulated data back here.
         self.CartCoord=self.manipulate.CartCoord
      
      #set threshold for smallest HR-factors:
      HRthresh=re.findall(r"(?<=HRthreshold=)[ \d.]+",self.opt,re.M)
      if HRthresh==[]:
         HRthresh=0.015
      else: 
         HRthresh=float(HRthresh[-1])

      #read HR-factors:
      self.readHR()
      
      #sort modes according to HR-factors
      index=np.argsort(self.HR, kind='heapsort')
      f=[]
      HR=[]
      for i in range(len(index)):
         HR.append(self.HR[index[i]])
         f.append(self.f[index[i]])
      self.f=[f,f]
      self.HR=HR

      if self.HR[-1]>=10: 
         #if larges HR-factor is too large
         loggingwrite(u'\nWARNING: THE HUANG-RHYS FACTOR SEEMS TO BE'+\
                                                         ' TOO LARGE !!\n')
         loggingwrite(u'        the spectrum will be calculated, but most probably'+\
                                       ' the input-state is inconsistent.\n')
      
      #print them
      self.log.write(u'HR-fact           freq     delta\n')
      for j in range(len(self.HR)-1,-1,-1):
         #select all 'big' HR-factors 
         if self.HR[j]>=HRthresh:
            self.log.write(u"%f   %f \n"%(self.HR[j], self.f[0][j]*self.Hartree2cm_1))

   def readHR(self):
      """ This function reads the HR-factors and electronic transition energy 
         from a given file and brings them into a 
         similar structure as they are used in the 'smallscript'.
         
         **PARAMETERS**
         logging:     This variable consists of two parts: logging[0] specifies 
                      the level of print-out (which is between 0- very detailed
                      and 4- only main information) and logging[1] is the file, 
                      already opened, to write the information in.
         hr_file:      the file where the information is found
         
         **RETURNS**
         initial:     a dummy-array that originally contains information about 
                      the inital states. Here at the moment only one
                      is allowed and only its length is relevant in the further 
                      programme.
         HRm:         a 2-dimensional array containing all Huang-Rhys-factors
         freqm:       analogously to HRm, containing the respective frequencies
         Energy:      the energy-difference of the electronic states. 
                      (in atomic units)
         
      """
      fi=open(self.hr_file, "r")
      f=mmap.mmap(fi.fileno(), 0, prot=mmap.PROT_READ)
      fi.close()
      Energy=re.findall(r"(?<=Delta E=)[ \d\.\-]*", f, re.I)
      assert len(Energy)==1, "Please specify only one energy-difference!"
      Energy=float(Energy[0])
      Energy/=self.Hartree2cm_1
      HR=[]
      funi=[]
      HRfreq=re.findall(r"HR-fact[\s]*freq[\s]*\n[\n\d\.\se\+\-]*", f, re.I)
      assert len(HRfreq)==1, "The file-format could not be read. exit now"
      HRf=re.findall(r"(?<=\n)[\d\.]*[\s]+[\d\.]*", HRfreq[0], re.I)
      for i in range(len(HRf)):
         line=re.findall(r"[\d.]+",HRf[i], re.I)
         HR.append(float(line[0]))
         funi.append(float(line[1])/self.Hartree2cm_1)
      #the following is just to be consistent with structure of 
      #                         HR calculated in first part
      self.f=funi
      self.HR=HR
      self.Energy[0]=Energy
      self.Energy[1]=0

   # function calcspect is enherited from FC_spects
   # more functions are not needed.

version='0.0.1'   
# End of HR_spects.py
