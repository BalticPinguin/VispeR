#!/usr/bin/python
#include [[functions_smsc.py]]
import functions_smsc as of 
import OPA 
import broaden as br
import sys, os, logging, re, mmap, numpy as np

def usage():
   print "usage: smallscript <input-file>"

def invokeLogging(mode="important"):
   if mode in ['all', 'ALL', 'All']:
      logging.basicConfig(filename='calculation.log',level=logging.DEBUG, format='%(message)s')
      logging.debug('use log-level all')
   elif mode in ['detailed', 'DETAILED', 'Detailed']:
      logging.basicConfig(filename='calculation.log',level=logging.INFO, format='%(message)s')
      logging.info('use log-level detailed')
   elif mode in ['medium', 'MEDIUM','Medium']:
      logging.basicConfig(filename='calculation.log',level=logging.WARNING, format='%(message)s')
      logging.warning('use log-level medium')
   elif mode in ['important', 'IMPORTANT','Important']:
      logging.basicConfig(filename='calculation.log',level=logging.ERROR, format='%(message)s')
   elif mode in ['short', 'SHORT','Short']:
      logging.basicConfig(filename='calculation.log',level=logging.CRITICAL, format='%(message)s')
   else:
      logging.basicConfig(filename='calculation.log',level=logging.ERROR, format='%(message)s')
      logging.error("logging-mode not recognized. Using 'important' instead")
   logging.info('Initializing the log-file')

def main(argv=None):
   assert len(argv)==1, 'exactly one argument required.'
   try:
      infile=open(argv[0], "r")
      f=mmap.mmap(infile.fileno(), 0, prot=mmap.PROT_READ)
      infile.close()
   except IOError:
      print "file", inputf, "not found."
      usage()
      return 2
   opts=[]
   todo=0
   if (re.search(r"HR-fact",f, re.I) is not None) is True:
      todo+=1
   opts.append(re.findall(r"(?<=HR-fact)[\w.,\(\) \=:\-]+", f, re.I))
   if (re.search(r"FC-spect",f, re.I) is not None) is True:
      if (re.search(r"HR-file: ",f, re.I) is not None) is True:
	 #calculation of HR-facts not neccecary
	 todo=2
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
   if todo>=16 or todo==0 or todo in [4,6,9,13]:
      print "options for calculation don't make sense. Please check the input-file!"
      return 0

   if np.mod(todo,2)==1: 
      #calculation up to HR-facts needed (FC- or Duschinsky spect)
      if opts[0]!=[]:
	 opt=opts[0][0]
      elif opts[1]!=[]:
	 opt=opts[1][0]
      elif opts[2]!=[]:
	 opt=opts[2][0]
      else:
	  print 'You want nothing to be calculated? Here it is:\n'
	  return 2
      #invoke logging: take options of HR-fact, if exists, else: from FC-spect, else: from Duschinsky-spect
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
	 invokeLogging()
      else:
	 invokeLogging(loglevel[0])
      initial=re.findall(r"(?<=initial: )[\w.]+",f, re.I)
      final=re.findall(r"(?<=final: )[\w.]+",f, re.I)
      HR, funi, Energy, J, K, f=of.CalculationHR(initial, final, opt)

   if np.mod(todo,4)>=2:
      #calculate FC-spect
      opt=opts[1][0]
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
	 invokeLogging()
      else:
	 invokeLogging(loglevel[0])
      if ((re.search(r"broaden",opt, re.I) is not None) is True) and todo<8:
	 todo+=8
      try: 
	 #test, whether HR-facts were calculated above
	 HR
      except NameError:
	 HRfile=re.findall(r"(?<=HR-file: )[\w.,\/\-]+",f, re.I)
	 assert len(HRfile)==1, 'There must be exactly one file specified containing HR-facts.'
	 HR, funi, E=of.ReadHR(HRfile[0])
	 Energy=np.zeros(2)
	 Energy[0]=E
	 Energy[1]=0
      T=re.findall(r"(?<=T=)[ \=\s\d\.]+", opt, re.M)
      if len(T)==0:
	 T=300
      else:
	 T=float(T[0])
      logging.info("temperature of system:"+repr(T))
      T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K
      part=re.findall(r"(particles:)[\s\d]*", opt, re.I)
      #############################################make this availible for multiple files!!
      if part==1:
	 linspect=of.calcspect(HR, funi, Energy[0]-Energy[1], 0, 5, 5, T, "OPA")
      elif part==2:
	 linspect=of.calcspect(HR, funi, Energy[0]-Energy[1], 0, 5, 5, T, "TPA")
      else:
	 linspect=of.calcspect(HR, funi, Energy[0]-Energy[1], 0, 5, 5, T)

   if np.mod(todo,8)>=4:
      #calculate Duschinsky-spect
      opt=opts[2][0]
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
	 invokeLogging()
      else:
	 invokeLogging(loglevel[0])
      if (re.search(r"broaden",opt, re.I) is not None) is True and todo<8: 
	 todo+=8
      try: #test, whether HR-facts were calculated above
	 J
      except NameError:
	 logging.critical('fatal error: No calculation of first part. But it is required')
	 return 2
      T=re.findall(r"(?<=T=)[ \=\s\d\.]+", opt, re.M)
      if len(T)==0:
	 T=300
      else:
	 T=float(T[0])
      logging.info("temperature of system:"+repr(T))
      T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K
      for i in range(len(initial)): #calculate separate line-spects for different states
	 linspect=OPA.FCfOPA(J[i], K[i], f, Energy[0]-Energy[1], 5, T, 0)

   if np.mod(todo,16)>=8:
      #calculate Broadening
      if opts[3]!=[]:
	 opt=opts[3][0]
      elif opts[2]!=[]:
	 opt=opts[2][0]
	 ############## take only sub-part of opt (where broaden) is specified...
	 ############## if not specified at all: take part of opts[1]
      elif opts[1]!=[]:
	 opt=opts[1][0]
	 ############## take only sub-part of opt (where broaden) is specified...
	 ############## if not specified at all: error!
      else:
	  print 'You want nothing to be calculated? Here it is:\n nothing'
	  return 2
      loglevel=re.findall(r"(?<=print\=)[\w]+",opt, re.I)
      if loglevel==[]:
	 invokeLogging()
      else:
	 invokeLogging(loglevel[0])
      T=re.findall(r"(?<=T=)[ \=\s\d\.]+", opt, re.M)
      if len(T)==0:
	 T=300
      else:
	 T=float(T[0])
      T*=8.6173324e-5/27.21138386 # multiplied by k_B in hartree/K

      try:
	 linspect ### does it have correct structure?
      except NameError:
	 linespectrum=re.findall(r"(?<=linspect: )[\w\.]+", f, re.I)
	 if linespectrum==[]:
	    linespectrum=re.findall(r"(?<=linespect:)[ \w\.]+", f, re.I)
	 assert len(linespectrum)==1, "if no spectrum calculation was done before, please specify a file containing line-spectrum."
	 freq=[]
	 intens=[]
	 mode=[]
	 with open(linespectrum[0]) as lines:
	    ##### handel: also files containing further information!
	    lis=[line.split() for line in lines]  # create a list of lists
	    for i,x in enumerate(lis):        # print the list items 
	       freq.append(float(x[0]))
	       intens.append(float(x[1]))
	       mode.append(float(x[2]))
	 linspect=np.zeros((3,len(freq)))
	 linspect[0]=np.matrix(freq)
	 linspect[1]=np.matrix(intens)
	 linspect[2]=np.matrix(mode)
      ################## change this to make it work with multiple files!!
      br.outspect(T, opt, linspect)
      ### if to nPA is specified: #### need energy-difference -> need to read it, if spectrum is taken from file...
	 #br.outspect(T, opt, linspect, Energy-difference)
   
if __name__ == "__main__":
   main(sys.argv[1:])
