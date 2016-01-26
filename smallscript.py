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

   #look, what kind of spectrum is to be obtained and initialise 
 
 # an object of respective class.
   if (re.search(r"model",f, re.I) is not None) is True:
      model = re.findall(r"model", f, re.I)[-1]
      if model == "FC":
         spect = FC-spect.init(f)
      if model == "CFC":
         spect = CFC-spect.init(f)
      if model == "GFC":
         spect = GFC-spect.init(f)
      if model == "HRs":
         spect = HR-spect.init(f)
      if model == "HRf":
         spect = HR-factors.init(f)
      else:
         spect = FC-spect.init(f)
   #INTRODUCTION END
   
    
if __name__ == "__main__":
   main(sys.argv[1:])

#version=1.6.1
# End of smallscript.py
