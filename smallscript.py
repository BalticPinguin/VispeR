#!/usr/bin/python2
# filename: smallscript.py

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
import time

def usage():
   print "usage: smallscript <input-file>"

def invokeLogging(logfile, mode="important"):
   """ initialises the logging-functionality
   ==PARAMETERS==
   logfile      name of file to be used as log-file. It is expected to be an array of 
                length 0 or one.
   mode:        5 different values are possible (see below); the lowest means: print 
                much, higher values mean 
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

def getTasks(f):
   """
   This function extracts the tasks together with all their options from the input-file.
   Additional text (that doesn't match the regexes) is simply ignored. 
   Therefore, be careful with typos.

   **PARAMETERS**
   f      the input-file (already opened)

   **RETURNS**
   opts   array containing all options of respective sub-tasks, which are 
          evaluated in the respective part
   todo   specifies the sub-tasks to be done in numeral values (powers of 2).
   """
   opts=[]
   todo=0
   # here: evaluate the file with respect to the tasks to be done
   # START PARSING TASKS
   if (re.search(r"HR-fact",f, re.I) is not None) is True:
      todo+=1
   opts.append(re.findall(r"(?<=HR-fact)[\w.,\(\) \=;:-]+", f, re.I))
   if (re.search(r"FC-spect",f, re.I) is not None) is True:
      if (re.search(r"HR-file: ",f, re.I) is not None) is True:
         #calculation of HR-facts not neccecary
         todo+=2
      else: #if 
         todo=3
   opts.append(re.findall(r"(?<=FC-spect)[\w\d\.\=,\(\):; -]+",f,re.I))
   if (re.search(r"Duschinsky-spect",f, re.I) is not None) is True:
      if todo==0:
         todo=5
      else:
         todo+=4
   opts.append(re.findall(r"(?<=Duschinsky-spect)[\w:\d\=.\(\),; -]+",f,re.I))
   if ((re.search(r"Broadening",f, re.I) is not None) is True) or\
       ((re.search(r"broaden",f, re.I) is not None) is True):
      todo+=8
   opts.append(re.findall(r"(?<=Broadening)[\w\d\.,:\(\)\=; -]+",f,re.I))
   # END PARSING TASKS
   
   # check, if the combination of input-arguments makes sense. This check is mainly for 
   # debugging purpose. There should be no invalid combination!
   if todo>=16 or todo in [0,4,6,9]: 
      print "options for calculation don't make sense. Please check the input-file!"
      print opts
   return opts, todo

def main(argv=None):
   """ This is the main-function of smallscript. 
       Its input-argument is the file containing all options. Here, it is evaluated
       wrt. the main tasks. 
   """
   #INTRODUCTION START
   assert len(argv)==1, 'exactly one argument required.'
   #open input-file (if existent and readable) and map it to f
   
   # try to open the input-file. If it doesn't exist or one is not alowed to open it,
   # print a usage-information and quit calculation.
   try:
      infile=open(argv[0], "r")
      f=mmap.mmap(infile.fileno(), 0, prot=mmap.PROT_READ)
      infile.close()
   except IOError:
      print "file", inputf, "not found."
      usage()
      return 2

   #If the input-file exists, get tasks and their options:
   opts, todo=getTasks(f)

   #invoke logging (write results into file specified by 'out: ' or into 'calculation.log')
   logfile=re.findall(r"(?<=out: )[\.\-_ \w]+",f, re.I)
   if logfile==[]:
      log=open("calculation.log", "a")
   else:
      log=open(logfile[0], "a")
   
   # now,  write the header to the output-file.
   log.write("\n==================================================================\n"
             "===================  output of smallscript  ======================\n"
             "==================================================================\n\n")
   log.write("   INPUT-FILE:\n")
   log.write(f)
   log.write(" \n   END OF INPUT-FILE \n\n")
   log.write("calculations to be done: %s\n"%(todo))
   log.close()
   #INTRODUCTION END
   
   # TASK EXECUTION START
   if np.mod(todo,2)==1: 
      #calculation of HR-facts needed (FC- or Duschinsky spect)
      #first find, where the options for this task is written
      if opts[0]!=[]:
         opt=opts[0][0]
      elif opts[1]!=[]:
         opt=opts[1][0]
      elif opts[2]!=[]:
         opt=opts[2][0]
      else: 
          # this is just error handling. Comes here only when input is dump
          print 'You want nothing to be calculated? Here it is:\n \n'
          return 2
      #invoke logging: take options of HR-fact, if exists, 
      #       else: from FC-spect, else: from Duschinsky-spect
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
         logging=invokeLogging(logfile)
      else:
         logging=invokeLogging(logfile,loglevel[0])
      # get file names
      final=re.findall(r"(?<=final: )[\w.\-]+",f, re.I)
      initial=re.findall(r"(?<=initial: )[\w.\-]+",f, re.I)
      HRthresh=re.findall(r"(?<=HRthreshold=)[ \d.]+",opt,re.I)
      if HRthresh==[]:
         # default threshold. Below this, normally no effects are visible.
         HRthresh=0.015
      else:
         HRthresh=float(HRthresh[0])
      # chose method that should be used for this
      method=re.findall(r"(?<=method:)[ \w]+",opt, re.I)
      if method==[]:
         logging[1].write("\nUse default method for the calculation of"+
                                                  " all quantities.\n")
         HR, funi, Energy, J, K, f=of.CalculationHR(logging, initial, 
                                                 final, opt, HRthresh)
      else:
         method=re.findall(r"[\w]+",method[0], re.I) # clean by spaces
         if method[0] in ["gradient", "Gradient", 'grad', 
                                    "gradient ", "grad "]:
            logging[1].write("\nUse method %s for the calculation of all quantities.\n"%(method[0]))
            HR, funi, Energy, J, K, f=of.gradientHR(logging, initial, 
                                                final, opt, HRthresh)
         elif (method[0] in ["shift", "SHIFT", "Shift"]) or\
              (method[0] in ["changed", "CHANGED", "Changed"]):
            logging[1].write("\nUse method %s for the calculation of all quantities.\n"%(method[0]))
            HR, funi, Energy, J, K, f=of.CalculationHR(logging, initial, 
                                                   final, opt, HRthresh)
         else:
            logging[1].write("method %s not recognised. Use Shift instead.\n"%(method))
            HR, funi, Energy, J, K, f=of.CalculationHR(logging, initial, 
                                                   final, opt, HRthresh)

   if np.mod(todo,4)>=2:
      #calculate FC-spect
      #here exists only one possibility for the options 
      opt=opts[1][0] 
      #set up logging
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
         try:
            logging[1].close()
            logging=invokeLogging(logfile)
         except UnboundLocalError:
            logging=invokeLogging(logfile)
      else:
         #i.e. if not specified: use default
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
      # get temperature of the system
      T=re.findall(r"(?<=T=)[ \=\s\d\.]+", opt, re.M)
      if len(T)==0:
         #use default
         T=300
      else:
         T=float(T[0])
      if logging[0]<=1:
         logging[1].write("temperature of system: "+repr(T)+"\n")
      T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K
      #The states of interest can be given in several forms:
      #   As 1 number (specifying the final-state vibrations that should be scanned
      #   As 2 numbers specifying the # of modes to be scanned in both states.
      states=re.findall(r"(?<=states=)[\d ,]*", opt, re.I)
      if len(states)==0:
         # i.e. no states given
	 states1=5
         states2=0
      else:
	 try:
            #i.e. 1 number given (else, the 'int' will give an error)
	    states1=int(states[0])
	    states2=0
	    logging[1].write("number of states: %d and %d \n"%(states1, states2))
	 except ValueError:
            try:
               # two numbers given, try to extract them
               states1=int(states[0].split(",")[0])
               states2=int(states[0].split(",")[1])
               logging[1].write("number of states: %d and %d \n"%(states1, states2))
            except ValueError:
               #unknown format. Use default and give a respective message.
               states1=5
               states2=0
               logging[1].write("!!number of vibrational states {0} is not an integer.",
                                    " Use default instead.\n".format(states1, states2))
    
      #Chose the method for calculating the FC-spectrum
      assert len(HR[0])>0, "There is no Huang-Rhys factor involved in the spectrum. You need at least one shifted mode!"
      if (re.search(r"changed", opt, re.I) is not None):
         logging[1].write("Calculate the stick-spectrum in FC-picture with %s excitations"\
            " and changing frequencies\n" %(states1) )
         linspect=of.changespect(logging, HR[0], funi, Energy[0]-Energy[1], 0, states1, states2, T)
      else:
         logging[1].write("Calculate the line-spectrum in FC-picture with %s and %s excitations.\n"%(states1,states2) )
         linspect=of.calcspect(logging, HR[0], funi[1], Energy[0]-Energy[1], 0, states1, states2, T)
      if ((re.search(r"broaden",opt, re.I) is not None) is True) and todo<8:
         if np.mod(todo,16)<8:
            #this is reached, if the option 'broaden' is specified in the FC-part.
            todo+=8

   if np.mod(todo,8)>=4:
      #calculate Duschinsky-spect
      opt=opts[2][0]

      #set the logging part here.
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
      
      # if the option 'broaden' is given in the Duschinsky-part, add it 
      #     to the todo-list
      if (re.search(r"broaden",opt, re.I) is not None) is True and todo<8: 
         if np.mod(todo,16)<8:
            todo+=8
     
      try: #test, whether HR-facts were calculated above
         J
      except NameError:
         # Here, it is not possible to read data from file.
         # And it is not worth implementing this because nearly no time is saved.
         logging[1].write('fatal error: No calculation of first part. But it is required')
         logging[1].close()
         return 2

      #set temperature of the system
      T=re.findall(r"(?<=T=)[ \=\s\d\.]+", opt, re.M)
      if len(T)==0:
         T=300
      else:
         T=float(T[0])
      if logging[0]<=1:
         logging[1].write("temperature of system: {0}\n".format(T))
      T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K
      states=re.findall(r"(?<=states=)[ \d]+", opt, re.I)
      if len(states)==0: # i.e. there was no specification by user
	 states=5
      else: 
         # try to read it
	 try:
	    states=int(states[0])
	    logging[1].write("number of states: {0}\n".format(states))
         # if unreadable: take default, but write some comment
	 except ValueError: 
	    logging[1].write("number of vibrational states {0} is not an integer. "
                              "Use default instead.\n".format(states))
	    states=5

      model=re.findall(r"(?<=model\=)[\w]+",opt, re.I)
      try: 
         model=model[0]
      except IndexError:
         #if not, than calculate in resort-model (default).
         for i in range(len(initial)): 
            k=[0,i]
            logging[1].write("\n Use the model 'resort' since none has been specified.\n")
            linspect=OPA.resortFCfOPA(logging, J, K, f[k], Energy[0]-Energy[1], states, T, 0)
      #else: try to find the model, the user wants to use
      if model in ['Simple', 'simple', 'SIMPLE']:
         logging[1].write("\n Use the model %s for the calculation of linespectrum.\n"%(model))
         linspect=OPA.simpleFCfOPA(logging, J, K, f, Energy[0]-Energy[1], states, T, 0)
      elif model in ['Resort', 'resort', 'RESORT']:
         logging[1].write("\n Use the model %s for the calculation of linespectrum.\n"%(model))
         linspect=OPA.resortFCfOPA(logging, J, K, f, Energy[0]-Energy[1], states, T, 0)
      elif model in ['Distributing', 'distributing', 'DISTRIBUTING', 'dist', 'DIST', 'Dist']:
         assert 1==2, "This model turned out to give really wrong results. It is even not good in any approximation Please use an other one."
         shifts=re.findall(r"(?<=maxshift\=)[\d]+",opt, re.I)
         if len(shifts)==1:
            logging[1].write("\n Use the model %s for the calculation of linespectrum.\n"
                              "Number of shifts taken into account: %d \n"%(model,int(shifts[0])))
            linspect=OPA.distFCfOPA(logging, J, K, f, Energy[0]-Energy[1], states, T, 0, int(shifts[0]))
         else:
            logging[1].write("\n Use the model %s for the calculation of linespectrum.\n"
                              "Number of shifts not specified, use 6 as default.\n"%(model) )
            linspect=OPA.distFCfOPA(logging, J, K, f, Energy[0]-Energy[1], states, T, 0, 6)
         # the threshold (4) can be made to be a parameter as well
      elif model in ["Unrestricted", 'UNRESTRITED', 'unrestricted', 'unrest']:
         #make 5 (number of excitations), 10 (number of vibrational mode taken into account) to parameters
         modes=re.findall(r"(?<=modes\=)[\d]+",opt, re.I)
         if len(modes)==1:
            logging[1].write("\n Use the model %s for the calculation of linespectrum.\n"
                  "Number of modes to be taken into account:  %s .\n"%(model,int(modes[0])))
            linspect=DR.unrestricted(logging, J, K, f, Energy[0]-Energy[1], states, T, 0, int(modes[0]))
         else:
            logging[1].write("\n Use the model %s for the calculation of linespectrum.\n"
                  "Number of modes not specified, use 10 as default.\n"%(model))
            linspect=DR.unrestricted(logging, J, K, f, Energy[0]-Energy[1], states, T,0, 10)
      else:
         # i.e.: the model was specified but not found. (give warning and do default)
         logging[1].write('An error occured. The option of "model" is not known! Please check the spelling,'\
               ' meanwile the Duschinsky-rotated spectrum is calculated using "resort".\n')
         linspect=OPA.resortFCfOPA(logging, J, K, f, Energy[0]-Energy[1], states, T, 0)

   
   np.set_printoptions(suppress=True)
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
      #test if there exists a line-spectrum already (calculated above
      try:
         linspect 
      except NameError:
         #if not, than extract it from the input-file specified by 'linspect: ' or 'linespect: '
         linespectrum=re.findall(r"(?<=linspect: )[\w\.]+", f, re.I)
         if linespectrum==[]:
            linespectrum=re.findall(r"(?<=linespect: )[\w\.]+", f, re.I)
         assert len(linespectrum)==1, "if no spectrum calculation was done before"+\
                                 ", please specify a file containing line-spectrum."
         freq=[]
         intens=[]
         mode=[]
         with open(linespectrum[0]) as lines:
            lis=[line.split() for line in lines]  # create a list of lists
            for i,x in enumerate(lis):            # print the list items 
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
      #this is real purpuse of this method; here the broadened spectrum is calculated.
      if re.search(r"broadenparallel", opt, re.I) is not None:
         br.parallelspect(logging, T, opt, linspect)
      else:
         print opt
         br.outspect(logging, T, opt, linspect)
      #if to nPA is specified: #### need energy-difference -> need to read it, if spectrum is taken from file...
      try:
         # if FC- and Dusch-spect were calculated; than probably both spectra need to be calculated in broadening...
         secondlinspect
         opt=opts[2][0]
         br.outspect(logging, T, opt, linspect)
      except NameError:
         opt=opts[0] #do something arbitrary
   # If the broadened spectrum is not calculated: print the line-spectrum
   # to not through away the results from that calculation. The line spect. will be 
   # written to the log-file than.
   else: 
      logging[1].write("intensity [arb. u.]    energy [cm^-1]\n")
      for i in range(len(linspect[0])):
         logging[1].write("%3.10f    \t%3.10f\n" %(linspect[1][i],linspect[0][i]))
    # TASK EXECUTION END
     
   # tell about the end of the programme
   logging[1].write("end of calculation reached. Normal exit.\n")
   logging[1].close()
    
if __name__ == "__main__":
   main(sys.argv[1:])

#version=1.6.1
# End of smallscript.py
