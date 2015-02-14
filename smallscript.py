#!/usr/bin/python
#include [[functions_smsc.py]]
import functions_smsc as of 
#include [[OPA.py]]
import OPA
#include [[Dusch_unrest.py]]
import Dusch_unrest as DR
#include [[broaden.py]]
import broaden as br
#include further dicts
import sys, re, mmap, numpy as np
from progressbar import ProgressBar
import time

def usage():
   print "usage: smallscript <input-file>"

def invokeLogging(logfile, mode="important"):
   """ initialises the logging-functionality
   ==PARAMETERS==
   logfile      name of file to be used as log-file. It is expected to be an array of length 0 or one.
   mode:        5 different values are possible (see below); the lowest means: print much, higher values mean 
                less printing
   ==RETURNS==
   logging:     number of respective mode
   log:         opened file to write in
   """
   if logfile==[]:
      log=open("calculation.log", "a")
   else:
      s=logfile[0].strip()
      log=open(s, "a")
   if mode in ['all', 'ALL', 'All', '0']:
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
      log.write("logging-mode not recognized. Using 'important' instead\n")
   return logging, log

def main(argv=None):
   assert len(argv)==1, 'exactly one argument required.'
   pbar = ProgressBar(maxval=100)
   #open input-file (if existent and readable) and map it to f
   try:
      infile=open(argv[0], "r")
      f=mmap.mmap(infile.fileno(), 0, prot=mmap.PROT_READ)
      infile.close()
   except IOError:
      print "file", inputf, "not found."
      usage()
      return 2
   #opts is an array containing all options of respective sub-tasks, which are evaluated in the respective part
   opts=[]
   # todo specifies the sub-tasks to be done in numerical values (powers of 2).
   todo=0
   # here: evaluate the file with respect to the tasks to be done
   if (re.search(r"HR-fact",f, re.I) is not None) is True:
      todo+=1
   opts.append(re.findall(r"(?<=HR-fact)[\w.,\(\) \=:\-]+", f, re.I))
   if (re.search(r"FC-spect",f, re.I) is not None) is True:
      if (re.search(r"HR-file: ",f, re.I) is not None) is True:
         #calculation of HR-facts not neccecary
         todo+=2
      else: #if 
         todo=3
   opts.append(re.findall(r"(?<=FC-spect)[\w\d\.\=,\(\): -]+",f,re.I))
   if (re.search(r"Duschinsky-spect",f, re.I) is not None) is True:
      if todo==0:
         todo=5
      else:
         todo+=4
   opts.append(re.findall(r"(?<=Duschinsky-spect)[\w:\d\=.\(\), -]+",f,re.I))
   if (re.search(r"Broadening",f, re.I) is not None) is True:
      todo+=8
   opts.append(re.findall(r"(?<=Broadening)[\w\d\.,:\(\)\= -]+",f,re.I))
   if todo>=16 or todo==0 or todo in [4,6,9]:
      print "options for calculation don't make sense. Please check the input-file!"
      print opts
      return 0
   logfile=re.findall(r"(?<=out: )[\.\-_ \w]+",f, re.I)
   if logfile==[]:
      log=open("calculation.log", "a")
   else:
      log=open(logfile[0], "a")
   log.write("calculations to be done: {0}\n".format(todo))
   log.close()

   pbar.update(5)
   if np.mod(todo,2)==1: 
      #calculation up to HR-facts needed (FC- or Duschinsky spect)
      #first find in where the options for this task are written in
      #print "start HR"
      if opts[0]!=[]:
         opt=opts[0][0]
      elif opts[1]!=[]:
         opt=opts[1][0]
      elif opts[2]!=[]:
         opt=opts[2][0]
      else:
          print 'You want nothing to be calculated? Here it is:\n \n'
          return 2
      #invoke logging: take options of HR-fact, if exists, else: from FC-spect, else: from Duschinsky-spect
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
         logging=invokeLogging(logfile)
      else:
         logging=invokeLogging(logfile,loglevel[0])
      initial=re.findall(r"(?<=initial: )[\w.]+",f, re.I)
      final=re.findall(r"(?<=final: )[\w.]+",f, re.I)
      ## the calculation of all needed quantities is done in this function
      method=re.findall(r"(?<=method: )[ \w]+",opt, re.I)
      if method==[]:
         HR, funi, Energy, J, K, f=of.CalculationHR(logging, initial, final, opt)
      elif method[0] in ["gradient", "Gradient", 'grad', "gradient ", "grad "]:
         ## test whether Duschinsky-rotation is needed
         print "gradient-method"
         HR, funi, Energy, K, f=of.gradientHR(logging, initial, final, opt)
      elif method[0] in ["shift", "SHIFT", "Shift"]:
         HR, funi, Energy, J, K, f=of.CalculationHR(logging, initial, final, opt)
      else:
         logging[1].write("method {0} not recognised. Use Shift instead.\n".format(method))
         HR, funi, Energy, J, K, f=of.CalculationHR(logging, initial, final, opt)
      for i in range(len(HR)):
         print "HR-fact,    frequency "
         for j in range(len(HR[0])):
            print HR[i][j], funi[i][j]*219474.63 

   pbar.update(12)
   if np.mod(todo,4)>=2:
      #calculate FC-spect
      #here exists only one possibility for the options 
      opt=opts[1][0] 
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
         try:
            logging[1].close()
            logging=invokeLogging(logfile)
         except UnboundLocalError:
            logging=invokeLogging(logfile)
      else:
         try:
            logging[1].close()
            logging=invokeLogging(logfile, loglevel[0])
         except UnboundLocalError:
            logging=invokeLogging(logfile, loglevel[0])
      try: 
         #test, whether HR-facts were calculated above
         HR
      except NameError:
         #otherwise they have to be extracted from file
         HRfile=re.findall(r"(?<=HR-file: )[\w.,\/\-]+",f, re.I)
         assert len(HRfile)==1, 'There must be exactly one file specified containing HR-facts.'
         initial, HR, funi, E=of.ReadHR(logging, HRfile[0])
         Energy=np.zeros(2)
         Energy[0]=E
         Energy[1]=0
      T=re.findall(r"(?<=T=)[ \=\s\d\.]+", opt, re.M)
      if len(T)==0:
         T=300
      else:
         T=float(T[0])
      if logging[0]<=1:
         logging[1].write("temperature of system: "+repr(T)+"\n")
      T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K
      states=re.findall(r"(?<=states=)[\d ]*", opt, re.I)
      if len(states)==0:
	 states=5
      else:
	 try:
	    states=int(states[0])
	    logging[1].write("number of states: {0}\n".format(states))
	 except ValueError:
	    logging[1].write("number of vibrational states {0} is not an integer. Use default instead.\n".format(states))
	    states=5
      for i in range(len(initial)):
         linspect=of.calcspect(logging, HR[i], funi[i], Energy[0]-Energy[1+i], 0, states, states, T)
      if ((re.search(r"broaden",opt, re.I) is not None) is True) and todo<8:
         if opts[2]!=[]:
            ## i.e.: the FC-spectrum has to be broadened and the Duschinsky-spect to be calculated
            secondlinspect=linspect
         if np.mod(todo,16)<8:
            todo+=8

   pbar.update(20)
   if np.mod(todo,8)>=4:
      #calculate Duschinsky-spect
      opt=opts[2][0]
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
         try:
            logging[1].close()
            logging=invokeLogging(logfile)
         except UnboundLocalError:
            logging=invokeLogging(logfile)
      else:
         try:
            logging[1].close()
            logging=invokeLogging(logfile, loglevel[0])
         except UnboundLocalError:
            logging=invokeLogging(logfile, loglevel[0])
      if (re.search(r"broaden",opt, re.I) is not None) is True and todo<8: 
         if np.mod(todo,16)<8:
            todo+=8
      try: #test, whether HR-facts were calculated above
         J
      except NameError:
         logging[1].write('fatal error: No calculation of first part. But it is required')
         logging[1].close()
         return 2
      T=re.findall(r"(?<=T=)[ \=\s\d\.]+", opt, re.M)
      if len(T)==0:
         T=300
      else:
         T=float(T[0])
      if logging[0]<=1:
         logging[1].write("temperature of system: {0}\n".format(T))
      T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K
      states=re.findall(r"(?<=states=)[ \d]+", opt, re.I)
      if len(states)==0:
	 states=5
      else:
	 try:
	    states=int(states[0])
	    logging[1].write("number of states: {0}\n".format(states))
	 except ValueError:
	    logging[1].write("number of vibrational states {0} is not an integer. Use default instead.\n".format(states))
	    states=5
      model=re.findall(r"(?<=model\=)[\w]+",opt, re.I)
      try:
         model=model[0]
      except IndexError:
         for i in range(len(initial)): 
            k=[0,i]
            linspect=OPA.resortFCfOPA(logging, J[i], K[i], f[k], Energy[0]-Energy[1], states, T, 0)
      if model in ['Simple', 'simple', 'SIMPLE']:
         for i in range(len(initial)):
            k=[0,i]
            ############# make linspect to append; at the moment it will be replaced!!!
            linspect=OPA.simpleFCfOPA(logging, J[i], K[i], f[k], Energy[0]-Energy[1], states, T, 0)
      elif model in ['Resort', 'resort', 'RESORT']:
         for i in range(len(initial)): #calculate separate line-spects for different states
            k=[0,i]
            linspect=OPA.resortFCfOPA(logging, J[i], K[i], f[k], Energy[0]-Energy[1], states, T, 0)
      elif model in ['Distributing', 'distributing', 'DISTRIBUTING', 'dist', 'DIST', 'Dist']:
         for i in range(len(initial)): #calculate separate line-spects for different states
            k=[0,i]
            shifts=re.findall(r"(?<=maxshift\=)[\d]+",opt, re.I)
            if len(shifts)==1:
               linspect=OPA.distFCfOPA(logging, J[i], K[i], f[k], Energy[0]-Energy[1], states, T, 0, int(shifts[0]))
            else:
               linspect=OPA.distFCfOPA(logging, J[i], K[i], f[k], Energy[0]-Energy[1], states, T, 0, 6)
            # the threshold (4) can be made to be a parameter as well
      elif model in ["Unrestricted", 'UNRESTRITED', 'unrestricted', 'unrest']:
         for i in range(len(initial)): #calculate separate line-spects for different states
            k=[0,i]
            #make 5 (number of excitations), 10 (number of vibrational mode taken into account) to parameters
            modes=re.findall(r"(?<=maxmodes\=)[\d]+",opt, re.I)
            if len(modes)==1:
               linspect=DR.unrestricted(logging, J[i], K[i], f[k], Energy[0]-Energy[1], states, T, 0, int(modes[0]))
            else:
               linspect=DR.unrestricted(logging, J[i], K[i], f[k], Energy[0]-Energy[1], states, T, 0, 10)
      else:
         logging[1].write('An error occured. The option of "model" is not known! Please check the spelling,'\
               ' meanwile the Duschinsky-rotated spectrum is calculated using "resort".\n')
         for i in range(len(initial)):
            k=[0,i]
            linspect=OPA.resortFCfOPA(logging, J[i], K[i], f[k], Energy[0]-Energy[1], states, T, 0)

   pbar.update(30)
   np.set_printoptions(suppress=True)
   #print 'linespect:', linspect.T
   if np.mod(todo,16)>=8:
      #calculate Broadening
      opt=0
      if opts[3]!=[]:
         opt=opts[3][0]
      else:
         for i in range(len(opts)):
            if opts[i]!=[]:
               if (re.search(r"(?<=broaden)[\w\.\-\= ,\(\):]", opts[i][0], re.M) is not None) is True:
                  opt=re.findall(r"(?<=broaden)[\w\.\-\= ,\(\):]+", opts[i][0], re.M)[0]
               break
      if opt==0:
         print 'You want nothing to be calculated? Here it is:\n nothing'
         return 2
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
         try:
            logging[1].close()
            logging=invokeLogging(logfile)
         except UnboundLocalError:
            logging=invokeLogging(logfile)
      else:
         try:
            logging[1].close()
            logging=invokeLogging(logfile, loglevel[0])
         except UnboundLocalError:
            logging=invokeLogging(logfile, loglevel[0])
      T=re.findall(r"(?<=T=)[ \=\s\d\.]+", opt, re.M)
      if len(T)==0:
         T=300
      else:
         T=float(T[0])
      T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K

      try:
         linspect 
      except NameError:
         linespectrum=re.findall(r"(?<=linspect: )[\w\.]+", f, re.I)
         if linespectrum==[]:
            linespectrum=re.findall(r"(?<=linespect: )[\w\.]+", f, re.I)
         assert len(linespectrum)==1, "if no spectrum calculation was done before"+\
                                 ", please specify a file containing line-spectrum."
         freq=[]
         intens=[]
         mode=[]
         with open(linespectrum[0]) as lines:
            ##### handel: also files containing further information!
            lis=[line.split() for line in lines]  # create a list of lists
            for i,x in enumerate(lis):        # print the list items 
               freq.append(float(x[0]))
               intens.append(float(x[1]))
               try:
                  mode.append(float(x[2]))
               except IndexError:
                  mode.append(42)
         linspect=np.zeros((3,len(freq)))
         linspect[0]=np.matrix(freq)
         linspect[1]=np.matrix(intens)
         linspect[2]=np.matrix(mode)
      ################## change this to make it work with multiple files!!
      pbar.update(50)
      br.outspect(logging, T, opt, linspect)
      pbar.update(80)
      ###if to nPA is specified: #### need energy-difference -> need to read it, if spectrum is taken from file...
      try:
         # if FC- and Dusch-spect were calculated; than probably both spectra need to be calculated in broadening...
         secondlinspect
         opt=opts[2][0]
         pbar.update(50)
         br.outspect(logging, T, opt, linspect)
         pbar.update(80)
      except NameError:
         opt=opts[0] #do something arbitrary

   logging[1].write("end of calculation reached. Normal exit.\n")
   logging[1].close()
   pbar.finish()
   
if __name__ == "__main__":
   main(sys.argv[1:])

version=1.3
